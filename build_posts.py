"""
Converts posts/*.md to posts/*.html using Pandoc and regenerates journal.json.
Run by GitHub Actions when posts/ changes; generated .html files are committed.
"""

import json
import subprocess
import datetime
import re
from pathlib import Path

import yaml

POSTS_DIR = Path("posts")
TEMPLATE = Path("_post-template.html")


def parse_frontmatter(text):
    m = re.match(r"^---\s*\n(.*?)\n---\s*\n", text, re.DOTALL)
    return yaml.safe_load(m.group(1)) if m else {}


posts = []

for md_file in sorted(POSTS_DIR.glob("*.md"), reverse=True):
    text = md_file.read_text(encoding="utf-8")
    meta = parse_frontmatter(text)

    slug = md_file.stem
    out_file = POSTS_DIR / f"{slug}.html"

    subprocess.run(
        [
            "pandoc", str(md_file),
            "--template", str(TEMPLATE),
            "--katex",
            "--highlight-style", "pygments",
            "-o", str(out_file),
        ],
        check=True,
    )
    print(f"Built {out_file}")

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
