#!/usr/bin/env python3
"""Build the static landing page into site/dist/.

Pure standard library. Reads:
  site/content.json  every claim on the page (edit this to change the text)
  CITATION.cff       author, dates, title, abstract, keywords, BibTeX
  site/style.css     inlined into the pages
  site/assets/       favicon.svg and og.png, copied as they are
  site/content.nightly.json  optional; merged over content.json for the nightly variant only
  site/nightly.css   styles of the nightly banner; never part of the default page

Usage:
  python site/build.py          build into site/dist/
  python site/build.py --pdf    also copy paper/ladder_short.pdf into site/dist/
                                (the Pages workflow does this copy itself)
  python site/build.py --variant nightly
                                build the work-in-progress preview into site/dist-nightly/.
                                It is served from the /nightly/ folder of the Pages site,
                                is noindex, has no citation tags and no JSON-LD, and
                                carries a banner with the git commit and the UTC build time.
  --out DIR                     write somewhere other than the default folder
  --commit HASH, --built-at ISO banner values of the nightly variant (default: this
                                checkout's HEAD and the current time; publish.py passes
                                them because it builds from an export that has no .git)

The build stops with a non-zero exit code if a sanity check fails: broken
JSON-LD, malformed XML or HTML, a dead internal link, a banned phrase, an
e-mail address, a private path, or a page heavier than the size budget. The
public page must carry the citation tags and no noindex; the nightly page the
reverse.
"""
import argparse
import datetime
import html
import json
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse

SITE = Path(__file__).resolve().parent
REPO = SITE.parent
DIST = SITE / "dist"
DIST_NIGHTLY = SITE / "dist-nightly"
NIGHTLY_PATH = "nightly/"  # where the preview lives, below the public page's address
BUDGET_BYTES = 150 * 1024  # page weight without the PDF and og.png

# Phrases the page must never carry. Checked case-insensitively against every
# generated text file. (A negation such as "not peer reviewed" is fine, so the
# review claim is checked as the pattern "peer-reviewed by" / "endorsed by".)
BANNED = [
    r"breakthrough",
    r"\bsolves?\b",
    r"proves quantum gravity",
    r"no equilibrium state exists for a ball",
    r"electric sectors balance",
    r"peer[- ]reviewed by",
    r"endorsed by",
    r"has been peer[- ]reviewed",
    r"\bis peer[- ]reviewed",
    r"\brefutes? (jacobson|senovilla|speranza|casini|longo)",
]
EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+\.[A-Za-z0-9.-]+")
PRIVATE_PATH = re.compile(r"[A-Za-z]:[\\/]+Users|/home/[a-z]|\\\\[a-z]+\\", re.I)


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------

def unquote(v):
    v = v.strip()
    if len(v) >= 2 and v[0] == v[-1] == '"':
        return v[1:-1]
    return v


def parse_cff(text):
    """Read the small subset of YAML that CITATION.cff uses."""
    data, cur = {}, None
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        stripped = raw.strip()
        if not raw.startswith(" ") and not stripped.startswith("-"):
            key, _, val = raw.partition(":")
            key, val = key.strip(), val.strip()
            if val == "":
                data[key], cur = [], key
            else:
                data[key], cur = unquote(val), None
        elif stripped.startswith("- "):
            item = stripped[2:]
            if cur == "authors" and ":" in item:
                k, _, v = item.partition(":")
                data[cur].append({k.strip(): unquote(v)})
            else:
                data[cur].append(unquote(item))
        elif cur == "authors" and ":" in stripped:
            k, _, v = stripped.partition(":")
            data[cur][-1][k.strip()] = unquote(v)
    for key in ("title", "abstract", "date-released", "version", "authors", "keywords", "license"):
        if key not in data:
            raise SystemExit(f"CITATION.cff has no '{key}'")
    return data


def fail(msg):
    raise SystemExit(f"build failed: {msg}")


def deep_merge(base, over, path="", added=None):
    """Overlay rule: objects merge key by key, arrays and scalars replace.

    Keys starting with "_" in the overlay are comments and are dropped. A key the
    overlay introduces that content.json lacks is recorded in `added`, because a
    typo would otherwise change nothing, silently.
    """
    out = dict(base)
    for k, v in over.items():
        if k.startswith("_"):
            continue
        here = f"{path}{k}"
        if k not in base:
            if added is not None:
                added.append(here)
            out[k] = v
        elif isinstance(base[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(base[k], v, here + ".", added)
        else:
            out[k] = v
    return out


def checkout_commit():
    """Short hash of this checkout's HEAD, with -dirty when tracked files differ."""
    try:
        h = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short=7", "HEAD"],
                           capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return h + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


# --------------------------------------------------------------------------
# Inline markup and block rendering
# --------------------------------------------------------------------------

class Renderer:
    LINK = re.compile(r"\[((?:[^\]`]|`[^`]*`)+)\]\(([^)\s]+)\)")

    def __init__(self, content, cff):
        self.c = content
        self.cff = cff
        s = content["site"]
        self.repo = s["repo"].rstrip("/")
        self.branch = s["branch"]
        self.macros = content["macros"]
        self.ids = set()

    # -- text ---------------------------------------------------------------
    def resolve(self, href):
        if href.startswith("blob:"):
            return f"{self.repo}/blob/{self.branch}/{href[5:]}"
        if href.startswith("tree:"):
            return f"{self.repo}/tree/{self.branch}/{href[5:]}"
        if href.startswith("arxiv:"):
            return f"https://arxiv.org/abs/{href[6:]}"
        return href

    def fmt(self, text):
        out = []
        for part in re.split(r"(`[^`]+`)", text):
            if len(part) > 1 and part.startswith("`") and part.endswith("`"):
                out.append(f"<code>{html.escape(part[1:-1], quote=False)}</code>")
                continue
            e = html.escape(part, quote=False)
            e = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", e)
            e = re.sub(r"\*(.+?)\*", r"<em>\1</em>", e)
            e = re.sub(r"(?<![\w])([A-Za-zℰℬ])_([A-Za-z]{1,2})(?![\w])", r"\1<sub>\2</sub>", e)  # c_l, x_i, ℰ_ij
            out.append(e)
        return "".join(out)

    def inline(self, text):
        out, pos = [], 0
        for m in self.LINK.finditer(text):
            out.append(self.fmt(text[pos:m.start()]))
            href = self.resolve(m.group(2))
            rel = ' rel="noopener"' if href.startswith("http") else ""
            out.append(f'<a href="{html.escape(href, quote=True)}"{rel}>{self.fmt(m.group(1))}</a>')
            pos = m.end()
        out.append(self.fmt(text[pos:]))
        return "".join(out)

    # -- math ---------------------------------------------------------------
    def mathml(self, src):
        def sub(m):
            if m.group(1) not in self.macros:
                fail(f"unknown math macro {m.group(1)}")
            return self.macros[m.group(1)]
        for _ in range(3):  # macros may use macros
            src = re.sub(r"\{\{(\w+)\}\}", sub, src)
        try:
            ET.fromstring(f"<math>{src}</math>")
        except ET.ParseError as e:
            fail(f"MathML is not well-formed ({e}): {src[:80]}")
        return src

    # -- blocks -------------------------------------------------------------
    def blocks(self, blocks):
        return "\n".join(self.block(b) for b in blocks)

    def block(self, b):
        t = b["type"]
        i = self.inline
        if t == "p":
            return f"<p>{i(b['text'])}</p>"
        if t == "h3":
            return f"<h3>{i(b['text'])}</h3>"
        if t in ("ul", "ol"):
            items = "".join(f"<li>{i(x)}</li>" for x in b["items"])
            return f"<{t}>{items}</{t}>"
        if t == "table":
            cols = b["columns"]
            head = "".join(f'<th scope="col">{i(c)}</th>' for c in cols)
            rows = []
            for r in b["rows"]:
                cells = [f'<th scope="row">{i(r[0])}</th>']
                for label, val in zip(cols[1:], r[1:]):
                    cells.append(f'<td data-label="{html.escape(label.replace("*", ""), quote=True)}"><span>{i(val)}</span></td>')
                rows.append("<tr>" + "".join(cells) + "</tr>")
            return (f'<div class="table-wrap"><table class="cards"><caption>{i(b["caption"])}</caption>'
                    f"<thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>")
        if t == "math":
            figs = []
            for it in b["items"]:
                figs.append(
                    '<figure class="eq">'
                    f'<figcaption>{i(it["label"])}</figcaption>'
                    f'<div class="math-scroll"><math display="block" alttext="{html.escape(it["alt"], quote=True)}">'
                    f'{self.mathml(it["mathml"])}</math></div></figure>')
            return "\n".join(figs)
        if t == "callout":
            hid = f"callout-{len(self.ids)}"
            self.ids.add(hid)
            return (f'<aside class="callout {b["kind"]}" aria-labelledby="{hid}">'
                    f'<h3 id="{hid}">{i(b["title"])}</h3>{self.blocks(b["blocks"])}</aside>')
        if t == "status":
            lis = []
            for it in b["items"]:
                text = it["text"]
                if text and text[0].isalpha():
                    text = text[0].upper() + text[1:]
                lis.append(f'<li><span class="tag {it["kind"]}">{i(it["label"])}</span>{i(text)}</li>')
            return f'<ul class="status">{"".join(lis)}</ul>'
        if t == "code":
            return f'<pre><code class="language-{html.escape(b["lang"])}">{html.escape(b["text"], quote=False)}</code></pre>'
        if t == "bibtex":
            return f'<pre><code class="language-bibtex">{html.escape(self.bibtex(b), quote=False)}</code></pre>'
        fail(f"unknown block type {t}")

    def bibtex(self, b):
        def esc(s):
            return re.sub(r"([&%$#_])", r"\\\1", s)
        authors = " and ".join(f"{a['family-names']}, {a['given-names']}" for a in self.cff["authors"])
        y, m, _ = self.cff["date-released"].split("-")
        month = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"][int(m) - 1]
        note = b["note"].replace("{version}", self.c["paper"]["version_number"])
        return (f"@misc{{{b['key']},\n"
                f"  author = {{{esc(authors)}}},\n"
                f"  title  = {{{{{esc(self.cff['title'])}}}}},\n"
                f"  year   = {{{y}}},\n"
                f"  month  = {month},\n"
                f"  note   = {{{esc(note)}}},\n"
                f"  url    = {{{self.c['site']['url']}}}\n"
                f"}}")


# --------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------

def attr(s):
    return html.escape(s, quote=True)


def json_for_script(obj):
    s = json.dumps(obj, ensure_ascii=False, indent=2)
    return s.replace("</", "<\\/").replace("<!--", "<\\!--")


def build_index(c, cff, css, r, nightly=None):
    """`nightly` is None for the public page, else {"public_url", "commit", "built"}."""
    s, p, h, hero = c["site"], c["paper"], c["head"], c["hero"]
    url = s["url"]
    author = f"{cff['authors'][0]['given-names']} {cff['authors'][0]['family-names']}"
    published = cff["date-released"]
    pdf_url = url + s["pdf"]
    og_url = url + "og.png"
    scholar_date = published.replace("-", "/")
    keywords = list(cff["keywords"]) + h["keywords_extra"]
    jl = c["jsonld"]

    jsonld = {
        "@context": "https://schema.org",
        "@graph": [
            {
                "@type": "ScholarlyArticle",
                "@id": url + "#paper",
                "name": cff["title"],
                "headline": hero["title"],
                "alternativeHeadline": hero["subtitle"],
                "abstract": cff["abstract"],
                "author": {"@type": "Person", "name": author},
                "datePublished": published,
                "dateModified": s["updated"],
                "inLanguage": s["language"],
                "url": url,
                "mainEntityOfPage": url,
                "image": og_url,
                "isAccessibleForFree": True,
                "license": jl["article_license"],
                "keywords": keywords,
                "version": p["version_number"],
                "creativeWorkStatus": jl["work_status"],
                "encoding": {
                    "@type": "MediaObject",
                    "contentUrl": pdf_url,
                    "encodingFormat": "application/pdf",
                },
                "isBasedOn": {"@id": url + "#code"},
            },
            {
                "@type": "SoftwareSourceCode",
                "@id": url + "#code",
                "name": jl["code_name"],
                "description": jl["code_description"],
                "codeRepository": s["repo"],
                "url": s["repo"],
                "programmingLanguage": "Python",
                "author": {"@type": "Person", "name": author},
                "datePublished": published,
                "dateModified": s["updated"],
                "version": cff["version"],
                "isAccessibleForFree": True,
                "license": jl["code_license"],
            },
        ],
    }

    metas = [
        ("name", "description", h["description"]),
        ("name", "robots", "index, follow, max-image-preview:large, max-snippet:-1"),
        ("name", "author", author),
        ("name", "citation_title", cff["title"]),
        ("name", "citation_author", author),
        ("name", "citation_publication_date", scholar_date),
        ("name", "citation_pdf_url", pdf_url),
        ("name", "citation_abstract_html_url", url),
        ("name", "citation_language", s["language"]),
        ("name", "citation_keywords", "; ".join(keywords)),
        ("property", "og:type", "website"),
        ("property", "og:site_name", s["name"]),
        ("property", "og:title", h["social_title"]),
        ("property", "og:description", h["social_description"]),
        ("property", "og:url", url),
        ("property", "og:image", og_url),
        ("property", "og:image:type", "image/png"),
        ("property", "og:image:width", "1200"),
        ("property", "og:image:height", "630"),
        ("property", "og:image:alt", h["og_alt"]),
        ("property", "og:locale", "en"),
        ("name", "twitter:card", "summary_large_image"),
        ("name", "twitter:title", h["social_title"]),
        ("name", "twitter:description", h["social_description"]),
        ("name", "twitter:image", og_url),
        ("name", "twitter:image:alt", h["og_alt"]),
    ]
    if nightly:
        # The preview must never compete with the public page: no indexing, no Scholar tags.
        metas = [(k, n, "noindex, nofollow" if n == "robots" else v)
                 for k, n, v in metas if not n.startswith("citation_")]
    meta_html = "\n".join(f'<meta {k}="{attr(n)}" content="{attr(v)}">' for k, n, v in metas)
    ld_html = "" if nightly else f'<script type="application/ld+json">\n{json_for_script(jsonld)}\n</script>\n'
    banner = ""
    if nightly:
        css = css + "\n" + (SITE / "nightly.css").read_text(encoding="utf-8").strip()
        banner = (
            '<div class="nightly" role="note" aria-label="Nightly preview"><p class="wrap">'
            '<strong>Nightly preview:</strong> work in progress, may be wrong. '
            f'<a href="{attr(nightly["public_url"])}">The public page is here.</a> '
            f'<span class="build">Build <code>{html.escape(nightly["commit"])}</code>, '
            f'{html.escape(nightly["built"])}.</span></p></div>\n')

    nav = "".join(f'<li><a href="#{n["id"]}">{html.escape(n["label"])}</a></li>' for n in c["nav"])
    buttons = "".join(
        f'<li><a{" class=primary" if b.get("primary") else ""} href="{attr(b["href"])}">{html.escape(b["label"])}</a></li>'
        for b in hero["buttons"])
    sections = "\n".join(
        f'<section id="{sec["id"]}" aria-labelledby="{sec["id"]}-h">'
        f'<h2 id="{sec["id"]}-h">{r.inline(sec["heading"])}</h2>\n{r.blocks(sec["blocks"])}</section>'
        for sec in c["sections"])
    foot_links = "".join(f'<li><a href="{attr(l["href"])}">{html.escape(l["label"])}</a></li>' for l in c["footer"]["links"])

    return f"""<!doctype html>
<html lang="{s['language']}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(h['title'])}</title>
<link rel="canonical" href="{attr(url)}">
{meta_html}
<meta name="color-scheme" content="light dark">
<meta name="theme-color" content="#fbfaf6" media="(prefers-color-scheme: light)">
<meta name="theme-color" content="#131417" media="(prefers-color-scheme: dark)">
<link rel="icon" href="favicon.svg" type="image/svg+xml">
{ld_html}<style>
{css}
</style>
</head>
<body>
<a class="skip" href="#main">Skip to the content</a>
{banner}<header class="top">
<div class="wrap">
<p class="kicker">{r.inline(hero['kicker'])}</p>
<h1>{r.inline(hero['title'])}</h1>
<p class="subtitle">{r.inline(hero['subtitle'])}</p>
<p class="meta">{html.escape(author)}. {html.escape(p['version_label'])}. {html.escape(p['status_label'])}</p>
<p class="lead">{r.inline(hero['lead'])}</p>
<aside class="notice" aria-labelledby="notice-h"><h2 id="notice-h">{r.inline(hero['notice_title'])}</h2><p>{r.inline(hero['notice'])}</p></aside>
<ul class="buttons">{buttons}</ul>
<nav class="toc" aria-label="On this page"><ul>{nav}</ul></nav>
</div>
</header>
<main id="main" class="wrap">
{sections}
</main>
<footer>
<div class="wrap">
<ul>{foot_links}</ul>
<p>{r.inline(c['footer']['text'])}</p>
<p>Page updated {html.escape(s['updated'])}.</p>
</div>
</footer>
</body>
</html>
"""


def build_404(c, css, r):
    s = c["site"]
    base = urlparse(s["url"]).path or "/"
    nf = c["not_found"]
    return f"""<!doctype html>
<html lang="{s['language']}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(nf['title'])} | {html.escape(s['name'])}</title>
<meta name="robots" content="noindex">
<meta name="color-scheme" content="light dark">
<link rel="icon" href="{attr(base)}favicon.svg" type="image/svg+xml">
<style>
{css}
</style>
</head>
<body>
<main class="wrap" style="padding-top:3rem">
<h1>{html.escape(nf['title'])}</h1>
<p>{r.inline(nf['text'])}</p>
<ul class="buttons">
<li><a class="primary" href="{attr(base)}">Go to the graviton ladder</a></li>
<li><a href="{attr(base)}{attr(s['pdf'])}">Paper (PDF)</a></li>
<li><a href="{attr(s['repo'])}">Repository</a></li>
</ul>
</main>
</body>
</html>
"""


def build_sitemap(c):
    s = c["site"]
    def entry(loc):
        return f"  <url>\n    <loc>{html.escape(loc)}</loc>\n    <lastmod>{s['updated']}</lastmod>\n  </url>\n"
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
            + entry(s["url"]) + entry(s["url"] + s["pdf"]) + "</urlset>\n")


def build_robots(c):
    return f"User-agent: *\nAllow: /\n\nSitemap: {c['site']['url']}sitemap.xml\n"


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------

VOID = {"meta", "link", "br", "hr", "img", "input", "source"}


class Balance(HTMLParser):
    """Fail on mis-nested or unclosed tags; collect ids, links and sources."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.ids, self.refs, self.errors = [], set(), [], []

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if "id" in d:
            if d["id"] in self.ids:
                self.errors.append(f"duplicate id {d['id']}")
            self.ids.add(d["id"])
        for k in ("href", "src"):
            if k in d and d[k]:
                self.refs.append(d[k])
        if tag not in VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"unexpected </{tag}> (open: {self.stack[-3:]})")
            if tag in self.stack:
                while self.stack and self.stack.pop() != tag:
                    pass
        else:
            self.stack.pop()


def check(c, dist, nightly=False):
    problems = []
    base = urlparse(c["site"]["url"]).path or "/"

    idx = (dist / "index.html").read_text(encoding="utf-8")
    m = re.search(r'<script type="application/ld\+json">\n(.*?)\n</script>', idx, re.S)
    has_noindex = re.search(r'<meta name="robots" content="[^"]*noindex', idx) is not None
    has_citation = "citation_" in idx
    canonical = re.search(r'<link rel="canonical" href="([^"]*)"', idx)
    if nightly:
        # The preview must stay out of search engines and Google Scholar.
        if not has_noindex:
            problems.append("nightly index.html lacks the noindex robots meta")
        if has_citation:
            problems.append("nightly index.html carries citation_ tags")
        if m or "application/ld+json" in idx:
            problems.append("nightly index.html carries JSON-LD")
        if 'class="nightly"' not in idx:
            problems.append("nightly index.html has no banner")
        if not canonical or canonical.group(1) != c["site"]["url"]:
            problems.append(f"nightly canonical is not {c['site']['url']}")
    else:
        if "noindex" in idx:
            problems.append("public index.html carries noindex")
        if not has_citation:
            problems.append("public index.html lacks the citation_ tags")
        if 'class="nightly"' in idx:
            problems.append("public index.html carries the nightly banner")

    # JSON-LD parses and has the fields search engines need (the preview has none).
    if nightly:
        pass
    elif not m:
        problems.append("no JSON-LD block")
    else:
        try:
            graph = json.loads(m.group(1))["@graph"]
            by_type = {n["@type"]: n for n in graph}
            need = {"ScholarlyArticle": ["name", "author", "datePublished", "license", "isAccessibleForFree", "encoding"],
                    "SoftwareSourceCode": ["name", "author", "datePublished", "license", "codeRepository", "isAccessibleForFree"]}
            for typ, keys in need.items():
                for k in keys:
                    if k not in by_type.get(typ, {}):
                        problems.append(f"JSON-LD {typ} lacks {k}")
        except (ValueError, KeyError) as e:
            problems.append(f"JSON-LD does not parse: {e}")

    # Well-formedness.
    if (dist / "sitemap.xml").exists():
        try:
            ET.parse(dist / "sitemap.xml")
        except ET.ParseError as e:
            problems.append(f"sitemap.xml: {e}")
    for name in [n for n in ("index.html", "404.html") if (dist / n).exists()]:
        bal = Balance()
        bal.feed((dist / name).read_text(encoding="utf-8"))
        problems += [f"{name}: {e}" for e in bal.errors]
        if bal.stack:
            problems.append(f"{name}: unclosed tags {bal.stack}")
        # Internal links must exist.
        for ref in bal.refs:
            if re.match(r"[a-z][a-z0-9+.-]*:", ref, re.I):
                continue
            path, _, frag = ref.partition("#")
            if path.startswith(base) and base != "/":
                path = path[len(base):]
            path = path.lstrip("/")
            if path and not (dist / path).exists():
                if not (path == "ladder_short.pdf"):  # copied by the workflow, or by --pdf
                    problems.append(f"{name}: dead internal link {ref}")
            if frag and not path and name == "index.html" and frag not in bal.ids:
                problems.append(f"{name}: dead anchor #{frag}")

    # Banned phrases, e-mail, private paths, in every text file.
    for f in sorted(dist.iterdir()):
        if f.suffix in (".html", ".txt", ".xml", ".svg"):
            text = f.read_text(encoding="utf-8")
            for pat in BANNED:
                if re.search(pat, text, re.I | re.S):
                    problems.append(f"{f.name}: banned phrase /{pat}/")
            if EMAIL.search(text):
                problems.append(f"{f.name}: looks like an e-mail address: {EMAIL.search(text).group(0)}")
            if PRIVATE_PATH.search(text):
                problems.append(f"{f.name}: private path: {PRIVATE_PATH.search(text).group(0)}")

    # Size budget.
    weight = sum(f.stat().st_size for f in dist.iterdir()
                 if f.is_file() and f.name not in ("og.png", "ladder_short.pdf"))
    if weight > BUDGET_BYTES:
        problems.append(f"page weight {weight} bytes exceeds the {BUDGET_BYTES} budget")
    return problems, weight


# --------------------------------------------------------------------------

def load_content(nightly):
    """content.json, with content.nightly.json merged over it for the nightly variant."""
    content = json.loads((SITE / "content.json").read_text(encoding="utf-8"))
    if not nightly:
        return content
    public_url = content["site"]["url"]
    overlay_file = SITE / "content.nightly.json"
    if overlay_file.exists():
        try:
            overlay = json.loads(overlay_file.read_text(encoding="utf-8"))
        except ValueError as e:
            fail(f"content.nightly.json does not parse: {e}")
        added = []
        content = deep_merge(content, overlay, added=added)
        print(f"overlay: merged {overlay_file.name} over content.json")
        for key in added:
            print(f"overlay: warning, '{key}' is not in content.json (a typo, or a deliberate addition?)")
    # The preview's address is not up to the overlay: every absolute URL points at /nightly/.
    content["site"]["url"] = public_url + NIGHTLY_PATH
    return content, public_url


def main(argv):
    ap = argparse.ArgumentParser(description="Build the landing page (see the module docstring).")
    ap.add_argument("--pdf", action="store_true", help="copy the paper's PDF next to the page")
    ap.add_argument("--variant", choices=("default", "nightly"), default="default")
    ap.add_argument("--out", type=Path, help="output folder (default: site/dist, or site/dist-nightly)")
    ap.add_argument("--commit", help="commit shown in the nightly banner (default: this checkout's HEAD)")
    ap.add_argument("--built-at", help="UTC build time for the banner, ISO 8601 such as 2026-11-15T03:00:00Z")
    args = ap.parse_args(argv)
    nightly = args.variant == "nightly"

    if nightly:
        content, public_url = load_content(True)
        commit = args.commit or checkout_commit()
        if not re.fullmatch(r"[0-9a-f]{7,40}(-dirty)?|unknown", commit):
            fail(f"--commit is not a git hash: {commit}")
        if args.built_at:
            try:
                when = datetime.datetime.fromisoformat(args.built_at.replace("Z", "+00:00")).astimezone(datetime.timezone.utc)
            except ValueError:
                fail(f"--built-at is not an ISO 8601 time: {args.built_at}")
        else:
            when = datetime.datetime.now(datetime.timezone.utc)
        banner = {"public_url": public_url, "commit": commit, "built": when.strftime("%Y-%m-%d %H:%M UTC")}
    else:
        content, banner = load_content(False), None
    out = (args.out or (DIST_NIGHTLY if nightly else DIST)).resolve()

    cff = parse_cff((REPO / "CITATION.cff").read_text(encoding="utf-8"))
    css = (SITE / "style.css").read_text(encoding="utf-8").strip()
    r = Renderer(content, cff)

    out.mkdir(parents=True, exist_ok=True)
    for old in out.iterdir():  # empty it, but keep the folder itself (a server may be using it)
        shutil.rmtree(old) if old.is_dir() else old.unlink()

    (out / "index.html").write_text(build_index(content, cff, css, r, banner), encoding="utf-8", newline="\n")
    if not nightly:
        # Only the public page has these; the preview must not advertise itself.
        (out / "404.html").write_text(build_404(content, css, r), encoding="utf-8", newline="\n")
        (out / "sitemap.xml").write_text(build_sitemap(content), encoding="utf-8", newline="\n")
        (out / "robots.txt").write_text(build_robots(content), encoding="utf-8", newline="\n")
        (out / ".nojekyll").write_text("", encoding="utf-8")
    shutil.copyfile(SITE / "assets" / "favicon.svg", out / "favicon.svg")
    shutil.copyfile(SITE / "assets" / "og.png", out / "og.png")
    if args.pdf:
        shutil.copyfile(REPO / content["site"]["pdf_source"], out / content["site"]["pdf"])

    problems, weight = check(content, out, nightly)
    for f in sorted(out.iterdir()):
        try:
            shown = f.relative_to(SITE)
        except ValueError:
            shown = f
        print(f"{f.stat().st_size:>8}  {shown}")
    print(f"page weight without PDF and og.png: {weight} bytes (budget {BUDGET_BYTES})")
    if problems:
        print("\nchecks failed:")
        for p in problems:
            print("  -", p)
        return 1
    print("checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
