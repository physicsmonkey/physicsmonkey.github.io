"""
Converts posts/*.md to posts/*.html using Pandoc and regenerates journal.json.
Also converts supplementary *.md files inside posts/*/ subdirectories.
Run by GitHub Actions when posts/ changes; generated .html files are committed.
"""

import json
import subprocess
import datetime
import re
from pathlib import Path

import yaml

POSTS_DIR = Path("posts")
TEMPLATE      = Path("_post-template.html")
SUPP_TEMPLATE = Path("_supp-template.html")

PANDOC_BASE = [
    "--katex",
    "--highlight-style", "pygments",
    "--lua-filter", "pandoc-link-targets.lua",
]


def parse_frontmatter(text):
    m = re.match(r"^---\s*\n(.*?)\n---\s*\n", text, re.DOTALL)
    return yaml.safe_load(m.group(1)) if m else {}


def pandoc(src, template, out):
    subprocess.run(
        ["pandoc", str(src), "--template", str(template), *PANDOC_BASE, "-o", str(out)],
        check=True,
    )
    print(f"Built {out}")


# ── Main posts (posts/*.md) ──────────────────────────────────────────────────

posts = []

for md_file in sorted(POSTS_DIR.glob("*.md"), reverse=True):
    text = md_file.read_text(encoding="utf-8")
    meta = parse_frontmatter(text)
    slug = md_file.stem
    pandoc(md_file, TEMPLATE, POSTS_DIR / f"{slug}.html")
    posts.append({
        "title":       str(meta.get("title", slug)),
        "date":        str(meta.get("date", "")),
        "slug":        slug,
        "url":         f"posts/{slug}.html",
        "tags":        meta.get("tags") or [],
        "description": str(meta.get("description", "")),
    })

out = {
    "updated": datetime.datetime.utcnow().isoformat() + "Z",
    "posts":   posts,
}
with open("journal.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print(f"Wrote journal.json with {len(posts)} post(s)")

# ── Supplementary files (posts/*/**/*.md) ────────────────────────────────────
# Converts reports and calculation readmes to HTML using the lightweight
# supplementary template. Raw data files (.csv, .py, .json, .png) are served
# as-is and need no conversion.

for md_file in sorted(POSTS_DIR.glob("*/**/*.md")):
    pandoc(md_file, SUPP_TEMPLATE, md_file.with_suffix(".html"))
