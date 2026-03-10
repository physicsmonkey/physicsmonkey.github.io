"""
Fetches Brian Swingle's papers from the arXiv API and writes papers.json.
Run by GitHub Actions on a schedule; papers.json is then served statically.
"""

import json
import datetime
import urllib.request
import xml.etree.ElementTree as ET

URL = (
    "https://export.arxiv.org/api/query"
    "?search_query=au:Swingle_Brian"
    "&sortBy=submittedDate&sortOrder=descending&max_results=150"
)

NS = {
    "atom":   "http://www.w3.org/2005/Atom",
    "arxiv":  "http://arxiv.org/schemas/atom",
}

def text(el):
    return " ".join((el.text or "").split()) if el is not None else ""

with urllib.request.urlopen(URL) as resp:
    root = ET.fromstring(resp.read())

papers = []
for entry in root.findall("atom:entry", NS):
    papers.append({
        "title":    text(entry.find("atom:title",   NS)),
        "abstract": text(entry.find("atom:summary", NS)),
        "published": text(entry.find("atom:published", NS)),
        "url":      text(entry.find("atom:id",      NS)),
        "authors":  [text(a) for a in entry.findall("atom:author/atom:name", NS)],
    })

out = {
    "updated": datetime.datetime.utcnow().isoformat() + "Z",
    "papers":  papers,
}

with open("papers.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

print(f"Wrote {len(papers)} papers to papers.json")
