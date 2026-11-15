#!/usr/bin/env python3
"""Build both Pages trees and prepare one commit for the gh-pages branch.

Pure standard library, plus git on PATH. Nothing is written inside this repository.

  python site/publish.py                         stable = main, nightly = this checkout's HEAD
  python site/publish.py --stable origin/main --nightly origin/nightly
  python site/publish.py --push                  also push the prepared commit

What it does:
  1. Exports each ref with `git archive` into a temporary folder, so uncommitted or
     untracked files can never leak into a publish, and builds it with that ref's own
     site/build.py (the stable page is whatever main builds; the preview is whatever
     the nightly ref builds).
  2. Lays out one tree:  /  the public page (default build of --stable) with the PDF,
     404.html, sitemap.xml (public page only), robots.txt, .nojekyll;
     /nightly/  the preview (nightly build of --nightly) with its own PDF copy.
  3. Checks the tree: every internal link and every own-site URL resolves inside it,
     no e-mail address, private path or banned phrase, noindex only on the preview and
     citation_ tags only on the public page, each page under the weight budget.
  4. Fetches gh-pages into a throwaway clone and makes ONE commit on top of its history,
     "Publish stable <hash> and nightly <hash>". Without --push it stops there and
     prints the exact `git push` command. Nothing changes when the same two commits
     are already published (--force to commit anyway).

A nightly ref with uncommitted changes (only checked when it is this checkout's HEAD)
is refused unless --allow-dirty is given. Even then only committed files are published.
Exit status: 0 done, 1 a build, check or git step failed, 2 refused (bad ref, dirty).
"""
import argparse
import datetime
import importlib.util
import json
import os
import posixpath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from html.parser import HTMLParser
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PAGES_BRANCH = "gh-pages"
NIGHTLY_DIR = "nightly"
EXPORT_PATHS = ["site", "CITATION.cff", "paper/ladder_short.pdf"]
TEXT_SUFFIXES = (".html", ".txt", ".xml", ".svg")

ROBOTS = """# Crawlers read only the robots.txt at the HOST root (https://ni1ra.github.io/robots.txt).
# A file inside a project site like this one is ignored, so it cannot protect anything.
# The preview in /nightly/ is kept out of search results by its noindex meta tag. It is
# deliberately not disallowed here either: a crawler that may not fetch a page never
# sees the noindex on it.
User-agent: *
Allow: /

Sitemap: {url}sitemap.xml
"""


class Stop(Exception):
    def __init__(self, msg, code=1):
        super().__init__(msg)
        self.code = code


def say(msg=""):
    print(msg, flush=True)


# --------------------------------------------------------------------------
# git
# --------------------------------------------------------------------------

def git(*args, cwd=REPO, check=True):
    p = subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True, encoding="utf-8")
    if check and p.returncode != 0:
        raise Stop(f"git {' '.join(args)} failed ({p.returncode}): {(p.stderr or p.stdout).strip()}")
    return p.stdout.strip()


def resolve_ref(ref):
    p = subprocess.run(["git", "rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}"],
                       cwd=str(REPO), capture_output=True, text=True)
    if p.returncode != 0:
        raise Stop(f"'{ref}' is not a commit in this repository (did you run `git fetch`?)", 2)
    full = p.stdout.strip()
    return full, git("rev-parse", "--short=7", full)


def export(commit, dest):
    """Unpack the files the site needs from one commit, straight from git's object store."""
    dest.mkdir(parents=True)
    # No line-ending conversion: a Windows autocrlf setting must not change published bytes.
    proc = subprocess.Popen(["git", "-c", "core.autocrlf=false", "-c", "core.eol=lf",
                             "archive", "--format=tar", commit, "--", *EXPORT_PATHS],
                            cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tar:
        try:
            tar.extractall(dest, filter="data")
        except TypeError:  # Python without extraction filters
            tar.extractall(dest)
    err = proc.stderr.read().decode("utf-8", "replace").strip()
    if proc.wait() != 0:
        raise Stop(f"git archive {commit[:7]} failed: {err}")


def build(src, nightly, commit_short, built_at):
    args = [sys.executable, "-X", "utf8", str(src / "site" / "build.py")]
    if nightly:
        args += ["--variant", "nightly", "--commit", commit_short, "--built-at", built_at]
    env = dict(os.environ, PYTHONUTF8="1")
    p = subprocess.run(args, cwd=str(src), capture_output=True, text=True, encoding="utf-8", env=env)
    label = "nightly" if nightly else "stable"
    if p.returncode != 0:
        raise Stop(f"{label} build failed ({p.returncode}):\n{p.stdout}{p.stderr}")
    out = src / "site" / ("dist-nightly" if nightly else "dist")
    if not (out / "index.html").is_file():
        raise Stop(f"{label} build produced no {out.name}/index.html (does that ref's build.py know --variant?)")
    last = [l for l in p.stdout.strip().splitlines() if l][-2:]
    say(f"  {label} build: {' | '.join(last)}")
    return out


# --------------------------------------------------------------------------
# tree
# --------------------------------------------------------------------------

def pdf_names(src):
    """(published name, source path inside the export) as the ref's content.json says."""
    try:
        site = json.loads((src / "site" / "content.json").read_text(encoding="utf-8"))["site"]
        return site["pdf"], site["pdf_source"]
    except (OSError, ValueError, KeyError):
        return "ladder_short.pdf", "paper/ladder_short.pdf"


def need(path):
    if not path.is_file():
        raise Stop(f"missing {path}")
    return path


def lay_out(tree, stable_dist, stable_src, nightly_dist, nightly_src, public_url):
    tree.mkdir(parents=True)
    for name in ("index.html", "404.html", "favicon.svg", "og.png", "sitemap.xml"):
        shutil.copyfile(need(stable_dist / name), tree / name)
    name, source = pdf_names(stable_src)
    shutil.copyfile(need(stable_src / source), tree / name)
    (tree / "robots.txt").write_text(ROBOTS.format(url=public_url), encoding="utf-8", newline="\n")
    (tree / ".nojekyll").write_text("", encoding="utf-8")

    night = tree / NIGHTLY_DIR
    night.mkdir()
    for name in ("index.html", "favicon.svg", "og.png"):
        shutil.copyfile(need(nightly_dist / name), night / name)
    name, source = pdf_names(nightly_src)
    shutil.copyfile(need(nightly_src / source), night / name)


class PageRefs(HTMLParser):
    """Collect ids, href/src values and own-site URLs in meta content."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids, self.refs, self.meta_urls = set(), [], []

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if d.get("id"):
            self.ids.add(d["id"])
        for k in ("href", "src"):
            if d.get(k):
                self.refs.append(d[k])
        if tag == "meta" and d.get("content", "").startswith("http"):
            self.meta_urls.append(d["content"])


def exists_exact(tree, rel):
    """File exists with this exact spelling (GitHub Pages is case sensitive, Windows is not)."""
    cur = tree
    parts = [p for p in rel.split("/") if p]
    for i, part in enumerate(parts):
        if not cur.is_dir() or part not in os.listdir(cur):
            return False
        cur = cur / part
    return cur.is_file()


def load_build_module():
    spec = importlib.util.spec_from_file_location("graviton_build", REPO / "site" / "build.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_tree(tree):
    """Return a list of problems found in the finished tree."""
    b = load_build_module()
    problems = []
    root_html = (tree / "index.html").read_text(encoding="utf-8")
    night_html = (tree / NIGHTLY_DIR / "index.html").read_text(encoding="utf-8")

    m = re.search(r'<link rel="canonical" href="([^"]*)"', root_html)
    if not m:
        return ["index.html has no canonical link"]
    public_url = m.group(1)
    night_url = public_url + NIGHTLY_DIR + "/"
    base_path = re.sub(r"^https?://[^/]+", "", public_url) or "/"

    # Which page is which.
    if "noindex" in root_html:
        problems.append("index.html (public) carries noindex")
    if "citation_" not in root_html or "application/ld+json" not in root_html:
        problems.append("index.html (public) lacks citation_ tags or JSON-LD")
    if not re.search(r'<meta name="robots" content="[^"]*noindex', night_html):
        problems.append("nightly/index.html lacks the noindex robots meta")
    if "citation_" in night_html or "application/ld+json" in night_html:
        problems.append("nightly/index.html carries citation_ tags or JSON-LD")
    c2 = re.search(r'<link rel="canonical" href="([^"]*)"', night_html)
    if not c2 or c2.group(1) != night_url:
        problems.append(f"nightly/index.html canonical is not {night_url}")
    if f'href="{public_url}"' not in night_html:
        problems.append("nightly/index.html banner does not link to the public page")
    for name in ("index.html", "404.html", "sitemap.xml"):
        if NIGHTLY_DIR in (tree / name).read_text(encoding="utf-8").lower():
            problems.append(f"{name} (public) mentions the nightly preview")

    # Links.
    def target(ref, page_dir):
        """Path inside the tree that `ref` points at, None if external, False if it escapes."""
        ref = ref.split("#", 1)[0].split("?", 1)[0]
        if re.match(r"[a-z][a-z0-9+.-]*:", ref, re.I):
            if not ref.startswith(public_url):
                return None
            path = ref[len(public_url):]
        elif ref.startswith("//"):
            return None
        elif ref.startswith("/"):
            if not ref.startswith(base_path):
                return False
            path = ref[len(base_path):]
        elif ref == "":
            path = page_dir
        else:
            joined = posixpath.normpath(posixpath.join(base_path + page_dir, ref))
            if not joined.startswith(base_path.rstrip("/")):
                return False
            path = joined[len(base_path):] + ("/" if ref.endswith("/") else "")
        if path == "" or path.endswith("/"):
            path += "index.html"
        return path

    pages = {"index.html": ("", root_html), f"{NIGHTLY_DIR}/index.html": (f"{NIGHTLY_DIR}/", night_html)}
    for name, (page_dir, text) in pages.items():
        refs = PageRefs()
        refs.feed(text)
        for ref in refs.refs + refs.meta_urls:
            if ref.startswith("#"):
                if len(ref) > 1 and ref[1:] not in refs.ids:
                    problems.append(f"{name}: dead anchor {ref}")
                continue
            t = target(ref, page_dir)
            if t is None:
                continue
            if t is False:
                problems.append(f"{name}: link leaves the site: {ref}")
            elif not exists_exact(tree, t):
                problems.append(f"{name}: dead internal link {ref} (no {t} in the tree)")
    for name in ("404.html",):
        refs = PageRefs()
        refs.feed((tree / name).read_text(encoding="utf-8"))
        for ref in refs.refs:
            t = target(ref, "")
            if t not in (None, False) and not exists_exact(tree, t):
                problems.append(f"{name}: dead internal link {ref}")

    # Content: banned phrases, e-mail addresses, private paths in every text file.
    for f in sorted(tree.rglob("*")):
        if f.is_file() and f.suffix in TEXT_SUFFIXES:
            rel = f.relative_to(tree).as_posix()
            text = f.read_text(encoding="utf-8")
            for pat in b.BANNED:
                if re.search(pat, text, re.I | re.S):
                    problems.append(f"{rel}: banned phrase /{pat}/")
            if b.EMAIL.search(text):
                problems.append(f"{rel}: looks like an e-mail address: {b.EMAIL.search(text).group(0)}")
            if b.PRIVATE_PATH.search(text):
                problems.append(f"{rel}: private path: {b.PRIVATE_PATH.search(text).group(0)}")

    # Weight of each page without the PDF and og.png, and the PDFs themselves.
    def weight(folder):
        return sum(f.stat().st_size for f in folder.iterdir()
                   if f.is_file() and f.name not in ("og.png",) and f.suffix != ".pdf")
    for label, folder in (("public page", tree), ("nightly page", tree / NIGHTLY_DIR)):
        w = weight(folder)
        say(f"  weight of the {label} without PDF and og.png: {w} bytes (budget {b.BUDGET_BYTES})")
        if w > b.BUDGET_BYTES:
            problems.append(f"{label} weighs {w} bytes, over the {b.BUDGET_BYTES} budget")
        for pdf in folder.glob("*.pdf"):
            if not pdf.read_bytes().startswith(b"%PDF"):
                problems.append(f"{pdf.relative_to(tree).as_posix()} is not a PDF")
    return problems


def list_tree(tree):
    say("Files in the Pages tree:")
    total = 0
    for f in sorted(p for p in tree.rglob("*") if p.is_file()):
        size = f.stat().st_size
        total += size
        say(f"  {size:>8}  {f.relative_to(tree).as_posix()}")
    say(f"  {total:>8}  total")


# --------------------------------------------------------------------------
# the gh-pages commit
# --------------------------------------------------------------------------

def prepare_commit(tree, work, remote, message, force):
    """Throwaway clone of gh-pages with the tree committed on top. Returns (clone, committed)."""
    url = git("remote", "get-url", remote)
    name, email = git("config", "user.name", check=False), git("config", "user.email", check=False)
    if not name or not email:
        raise Stop("git user.name / user.email are not set in this repository; set them before publishing", 2)

    clone = work / "pages"
    has_branch = bool(git("ls-remote", "--heads", url, PAGES_BRANCH))
    if has_branch:
        git("clone", "--quiet", "--no-checkout", "--single-branch", "--branch", PAGES_BRANCH, url, str(clone), cwd=work)
    else:
        clone.mkdir()
        git("init", "--quiet", "-b", PAGES_BRANCH, cwd=clone)
        git("remote", "add", "origin", url, cwd=clone)
        say(f"  note: {remote} has no {PAGES_BRANCH} branch yet, so this commit starts it")
    git("config", "core.autocrlf", "false", cwd=clone)
    git("config", "user.name", name, cwd=clone)
    git("config", "user.email", email, cwd=clone)

    parent = git("rev-parse", "--short=7", "HEAD", cwd=clone) if has_branch else None
    if has_branch and not force and git("log", "-1", "--format=%s", cwd=clone) == message:
        say(f"  {PAGES_BRANCH} already carries '{message}'; nothing to publish (--force overrides)")
        return clone, False

    # The tree replaces the branch's files; history stays.
    for item in tree.iterdir():
        if item.is_dir():
            shutil.copytree(item, clone / item.name)
        else:
            shutil.copyfile(item, clone / item.name)
    git("add", "-A", cwd=clone)
    if has_branch and subprocess.run(["git", "diff", "--cached", "--quiet", "HEAD"], cwd=str(clone)).returncode == 0:
        say(f"  {PAGES_BRANCH} already has exactly these files; nothing to publish")
        return clone, False
    git("commit", "--quiet", "-m", message, cwd=clone)
    say(f"  committed {git('rev-parse', '--short=7', 'HEAD', cwd=clone)} on top of "
        f"{parent or 'nothing (first commit)'}: {message}")
    say(git("show", "--stat", "--format=", "HEAD", cwd=clone))
    return clone, True


# --------------------------------------------------------------------------

def main(argv):
    ap = argparse.ArgumentParser(description="Prepare the gh-pages commit (see the module docstring).")
    ap.add_argument("--stable", default="main", help="ref published at the site root (default: main)")
    ap.add_argument("--nightly", default="HEAD", help="ref published under /nightly/ (default: this checkout's HEAD)")
    ap.add_argument("--allow-dirty", action="store_true", help="publish although the checkout has uncommitted changes")
    ap.add_argument("--push", action="store_true", help="push the prepared commit to gh-pages")
    ap.add_argument("--force", action="store_true", help="commit even if gh-pages already carries these two commits")
    ap.add_argument("--remote", default="origin")
    ap.add_argument("--keep", action="store_true", help="keep the temporary folder after a push or a failure")
    args = ap.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    work = Path(tempfile.mkdtemp(prefix="graviton-pages-"))
    keep = args.keep
    try:
        stable, stable_short = resolve_ref(args.stable)
        nightly, nightly_short = resolve_ref(args.nightly)
        head = git("rev-parse", "HEAD")
        dirty = git("status", "--porcelain", "--untracked-files=normal")
        if nightly == head and dirty:
            if not args.allow_dirty:
                raise Stop("the nightly ref is this checkout's HEAD and the checkout has uncommitted changes:\n"
                           + dirty + "\nCommit them, or pass --allow-dirty (only committed files are published).", 2)
            say("warning: uncommitted changes are NOT part of the preview; only the commit is.")
        say(f"stable  {args.stable} -> {stable_short}")
        say(f"nightly {args.nightly} -> {nightly_short}")
        say(f"working folder: {work}")

        built_at = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        export(stable, work / "stable-src")
        export(nightly, work / "nightly-src")
        stable_dist = build(work / "stable-src", False, stable_short, built_at)
        nightly_dist = build(work / "nightly-src", True, nightly_short, built_at)

        content = json.loads((work / "stable-src" / "site" / "content.json").read_text(encoding="utf-8"))
        tree = work / "tree"
        lay_out(tree, stable_dist, work / "stable-src", nightly_dist, work / "nightly-src", content["site"]["url"])
        say()
        list_tree(tree)
        say()
        say("Checks on the tree:")
        problems = check_tree(tree)
        if problems:
            say("checks FAILED:")
            for p in problems:
                say(f"  - {p}")
            return 1
        say("  all checks passed")
        shutil.rmtree(work / "stable-src")
        shutil.rmtree(work / "nightly-src")

        say()
        message = f"Publish stable {stable_short} and nightly {nightly_short}"
        clone, committed = prepare_commit(tree, work, args.remote, message, args.force)
        if not committed:
            keep = args.keep
            return 0
        push = f'git -C "{clone}" push origin HEAD:{PAGES_BRANCH}'
        say()
        if args.push:
            say(f"pushing: {push}")
            say(git("push", "origin", f"HEAD:{PAGES_BRANCH}", cwd=clone))
            say("pushed. GitHub Pages usually serves it within a minute or two.")
        else:
            keep = True
            say("Stopped before pushing. To publish exactly this commit, run:")
            say()
            say(f"  {push}")
            say()
            say(f"(or run this script again with --push; delete {work} when done)")
        return 0
    except Stop as e:
        print(f"publish stopped: {e}", file=sys.stderr, flush=True)
        keep = keep or args.keep
        return e.code
    finally:
        if not keep:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
