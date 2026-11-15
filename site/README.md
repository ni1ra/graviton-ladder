# Landing page

A one-page static site for the repository, built for search engines and Google Scholar.
It is served by GitHub Pages from the `gh-pages` branch at two addresses:

| Address | What it shows | Search engines |
|---|---|---|
| <https://ni1ra.github.io/graviton-ladder/> | The public page: the default build of `main`. | Indexed, with Google Scholar tags and JSON-LD. |
| <https://ni1ra.github.io/graviton-ladder/nightly/> | The nightly preview: the newest work in progress, built from another branch. Open it on a phone to test before it becomes the front page. | `noindex, nofollow`, no Scholar tags, no JSON-LD. |

The preview carries a thin banner at the top ("Nightly preview: work in progress, may be wrong.
The public page is here.") with the git commit and the UTC time it was built from. Otherwise it is
the same page. Its canonical URL and every other absolute URL point at `/nightly/`.

## Rebuild locally

```bash
python site/build.py --pdf
```

Python 3 standard library only, nothing to install. The page is written to `site/dist/`;
open `site/dist/index.html` in a browser. `--pdf` also copies `paper/ladder_short.pdf` next
to it (`publish.py` does that copy itself, so the PDF is never committed twice).
The build stops with a non-zero exit code if JSON-LD, the sitemap or the HTML is broken, an
internal link is dead, a banned phrase or an e-mail address appears, or the page is over 150 KB
(PDF and `og.png` excluded). It also stops if the public page lacks the `citation_` tags or
carries `noindex`, or if the nightly page lacks `noindex` or carries `citation_` tags or JSON-LD.

## Preview locally

```bash
python site/build.py --pdf                       # public page into site/dist/
python site/build.py --variant nightly --pdf     # preview into site/dist-nightly/
python -m http.server 8000 --directory site/dist-nightly
```

Then open <http://localhost:8000/>. The nightly build shows this checkout's HEAD in the banner
(with `-dirty` when tracked files differ) and the current time. Both output folders are
git-ignored. `--out DIR` writes somewhere else.

## Nightly content overlay

`site/content.nightly.json`, if present, is merged over `content.json` for the nightly build only.
The rule: **objects merge key by key; arrays and scalars replace** (a replaced array is replaced
as a whole, not element by element). Keys starting with an underscore are comments and are
dropped. `site.url` cannot be overridden: the preview's address is fixed. A key that the overlay
adds but `content.json` lacks is printed as a warning, because a typo would otherwise change
nothing.

The overlay shipped today only changes `paper.version_label`, so the mechanism is exercised. Put
text that is not yet public here; when it is ready, move it into `content.json` and delete it from
the overlay. The overlay goes through the same checks as the rest (banned phrases, e-mail
addresses).

## Publish

```bash
python site/publish.py                                   # stable = main, nightly = this checkout's HEAD
python site/publish.py --stable origin/main --nightly origin/nightly
python site/publish.py --push                            # also push
```

Python standard library and `git` only. `publish.py`:

1. Exports each ref with `git archive` into a temporary folder outside the repository (so
   uncommitted and untracked files can never leak) and builds it with that ref's own `build.py`.
   A nightly ref that is this checkout's HEAD and has uncommitted changes is refused unless
   `--allow-dirty` is given; even then only the committed files are published.
2. Lays out the tree: the public page at the root with `ladder_short.pdf`, `404.html`,
   `sitemap.xml` (public page only), `robots.txt` and `.nojekyll`; the preview in `nightly/` with
   its own PDF copy.
3. Checks the tree: every internal link and own-site URL resolves; no e-mail address, private
   path or banned phrase; `noindex` only on the preview and `citation_` tags only on the public
   page; each page under 150 KB without PDF and `og.png`.
4. Fetches `gh-pages` into a throwaway clone and makes **one** commit on top of its history,
   `Publish stable <hash> and nightly <hash>`, then **stops** and prints the exact `git push`
   command. `--push` runs it. If the same two commits are already published it does nothing
   (`--force` commits anyway).

`robots.txt` is a courtesy only: crawlers read just the `robots.txt` at the host root
(`ni1ra.github.io/robots.txt`), never one inside a project site. The `noindex` meta tag is what
keeps the preview out of search results, and it is why `robots.txt` does not disallow `/nightly/`
(a crawler that may not fetch a page never sees its `noindex`).

## Keep the preview fresh (Windows scheduled task)

Actions cannot publish this (the login has no `workflow` scope), so a local task does.
`site/nightly.ps1` runs `git fetch`, then `publish.py --push` with `origin/main` and
`origin/nightly`, appends to `C:\Users\nira\work\artifacts\graviton-site\nightly-log.txt` (outside
the repository) and exits non-zero if either step fails. The preview is built from the branch
**`nightly`** on `origin`: push work there to see it on the phone after the next run. Without a
`nightly` branch the run fails and says so in the log.

Start it (every day at 03:00, as the logged-in user, while the PC is on):

```
schtasks /Create /TN "GravitonLadderNightly" /SC DAILY /ST 03:00 /TR "powershell.exe -NoProfile -ExecutionPolicy Bypass -File C:\Users\nira\work\graviton-ladder\site\nightly.ps1" /F
```

Stop it for good:

```
schtasks /Delete /TN "GravitonLadderNightly" /F
```

Run it once now: `schtasks /Run /TN "GravitonLadderNightly"`. Check it:
`schtasks /Query /TN "GravitonLadderNightly" /V /FO LIST`. A dry run that pushes nothing:
`powershell -File site\nightly.ps1 -DryRun`. Each run ends by itself within a minute or so, and
only ever pushes to `gh-pages`, never to `main`.

## Update the claims

Edit `site/content.json` only. Every sentence on the page, the page title and description, the
status list, the links and the paper's version label live there. Inline markup in strings:
`*emphasis*`, `**strong**`, `` `code` ``, `[text](url)`; `blob:path` and `tree:path` link into the
repository, `arxiv:ID` into arXiv. Formulas are MathML; `{{name}}` expands a macro from the
`macros` object.

Author, release date, title, abstract and keywords (page metadata, JSON-LD, the Google Scholar
tags and the BibTeX entry) are read from `CITATION.cff`, so change them there. When the paper
moves to a new version, update `paper.version_label`, `paper.version_number` and `site.updated`
in `content.json`, and the date and title in `CITATION.cff`.

## Other files

| File | What it is |
|---|---|
| `build.py` | Renders `content.json` into `dist/` (`--variant nightly`: into `dist-nightly/`). |
| `style.css` | The stylesheet, inlined into the pages. Light and dark follow the reader's system setting. |
| `nightly.css` | The banner's styles, inlined into the nightly page only. |
| `content.nightly.json` | Optional overlay for the nightly page (see above). |
| `publish.py` | Builds both pages and prepares the `gh-pages` commit. |
| `nightly.ps1` | The scheduled task's script. |
| `assets/favicon.svg` | Favicon. |
| `assets/og.png` | The 1200 x 630 social card, committed. Regenerate with `python site/make_og.py` (needs Pillow; the build does not). |

## Pages setup

Settings > Pages > Build and deployment > Source: **Deploy from a branch**, branch `gh-pages`,
folder `/ (root)`. Nothing builds on GitHub: the branch holds the finished files, and only
`publish.py` writes to it, one commit per publish, history kept. The `.nojekyll` file stops
GitHub running Jekyll over it.
