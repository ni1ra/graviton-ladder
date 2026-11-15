# Landing page

A one-page static site for the repository, built for search engines and Google Scholar.
It is served by GitHub Pages at <https://ni1ra.github.io/graviton-ladder/>.

## Rebuild locally

```bash
python site/build.py --pdf
```

Python 3 standard library only, nothing to install. The page is written to `site/dist/`;
open `site/dist/index.html` in a browser. `--pdf` also copies `paper/ladder_short.pdf` next
to it (the deploy workflow does that copy itself, so the PDF is never committed twice).
The build stops with a non-zero exit code if JSON-LD, the sitemap or the HTML is broken, an
internal link is dead, a banned phrase or an e-mail address appears, or the page is over 150 KB
(PDF and `og.png` excluded).

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
| `build.py` | Renders `content.json` into `dist/`. |
| `style.css` | The only stylesheet, inlined into the pages. Light and dark follow the reader's system setting. |
| `assets/favicon.svg` | Favicon. |
| `assets/og.png` | The 1200 x 630 social card, committed. Regenerate with `python site/make_og.py` (needs Pillow; the build does not). |

## Turning Pages on (once)

Settings > Pages > Build and deployment > Source: **GitHub Actions**. Then run the `pages`
workflow once from the Actions tab, or push a change under `site/`. The workflow then runs on every push to `main` that touches `site/`,
`paper/`, `CITATION.cff` or the workflow.
