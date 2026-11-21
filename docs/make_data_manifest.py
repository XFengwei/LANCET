"""Write data/MANIFEST.tsv (run from the repository root): every delivered file,
the top two levels of data/reduction and data/extraction, with origin,
reproducibility and the tasks that read each path."""
import os, re, subprocess
from pathlib import Path

root = Path.cwd()
readers = {}
texts = {}
for f in list(root.glob("tasks/*/task.toml")) + list(root.glob("tasks/*/*.py")) + list(root.glob("design/*/*.py")):
    if "paper_ii_referee_response" in str(f):
        continue
    owner = f.parent.name
    texts.setdefault(owner, "")
    texts[owner] += f.read_text()

def used_by(rel):
    names = sorted(o for o, t in texts.items() if rel in t)
    if not names:  # folder-level declarations (e.g. VVV tile directory, getsf folders)
        parent = rel
        while "/" in parent:
            parent = parent.rsplit("/", 1)[0]
            hits = sorted(o for o, t in texts.items() if f'"{parent}"' in t)
            if hits:
                return ",".join(hits)
    return ",".join(names)

def classify(rel):
    p = rel.removeprefix("data/delivered/")
    if p.startswith("alma/continuum/allscale-archived/"):
        return "derived", "no", "superseded v1-v3 12m+7m+TP combinations (Paper II); replaced by combine_continuum"
    if p.startswith("alma/continuum/IRDC316.8.TP7MTM2TM1.") or p.startswith("alma/continuum/IRDC316.8.TP7M."):
        if "TP7MTM2TM1" in p:
            return "derived", "task:combine_continuum", "promoted J-comb product"
        return "derived", "no", "older TP+7m combination"
    if p.startswith("alma/continuum/casa-archived/"):
        return "reduced", "reduction", "CASA 12m+7m, archived, not for science"
    if p.startswith("alma/"):
        return "reduced", "reduction", ""
    if p.startswith("archive/planck/") and not p.startswith("archive/planck/hfi_cutouts") and not p.startswith("archive/planck/PLCK"):
        return "derived", "task:planck_polarization", "earlier polarization products; planck_polarization rebuilds them"
    if p.startswith("archive/radio/subcube_G316.hdr"):
        return "derived", "no", "cut-out header (provenance)"
    if p.startswith("archive/paper-1/"):
        return "collaborator", "no", "published Paper I table"
    if p.startswith("archive/"):
        return "archive", "download", ""
    if p.startswith("sedfitter/"):
        return "derived", "no", "Paper I SED fitting (outside this repository)"
    if p.startswith("simulations/"):
        return "collaborator", "no", ""
    if p.startswith("total_power/"):
        return ("archive", "download", "SEDIGISM") if "SEDIGISM" in p else ("observed", "no", "")
    if p == "README.md" or p.endswith("README.md"):
        return "documentation", "", ""
    return "unknown", "", "check"

def size(path):
    if path.is_file():
        return path.stat().st_size
    out = subprocess.run(["du", "-sk", str(path)], capture_output=True, text=True).stdout.split()[0]
    return int(out) * 1024

rows = []
for path in sorted((root / "data/delivered").rglob("*")):
    if path.is_file() and path.name != ".DS_Store":
        rel = str(path.relative_to(root))
        origin, repro, note = classify(rel)
        if path.name.endswith("README.md"):
            origin, repro, note = "documentation", "", ""
        rows.append((rel, size(path), origin, repro, used_by(rel), note))
for base, origin, repro, note in (("data/reduction", "reduced", "reduction", "workspace (CASA, MIRIAD); not an analysis input unless declared"),
                                   ("data/extraction", "derived", "no", "getsf and astrodendro extractions; getsf is not rerun (D7)")):
    for level1 in sorted((root / base).iterdir()):
        if level1.name == ".DS_Store":
            continue
        children = [c for c in sorted(level1.iterdir()) if c.name != ".DS_Store"] if level1.is_dir() else []
        for item in ([level1] if not children else children):
            rel = str(item.relative_to(root)) + ("/" if item.is_dir() else "")
            rows.append((rel, size(item), origin if not item.name.endswith("README.md") else "documentation",
                         repro if not item.name.endswith("README.md") else "", used_by(rel.rstrip("/")), note if not item.name.endswith("README.md") else ""))
header = "path\tsize_bytes\torigin\treproducible\tused_by\tnote\n"
(root / "data/MANIFEST.tsv").write_text(header + "".join("\t".join(map(str, r)) + "\n" for r in rows))
from collections import Counter
print(len(rows), "rows;", Counter(r[2] for r in rows))
print("unknown:", [r[0] for r in rows if r[2] == "unknown"])
print("delivered files read by a task:", sum(1 for r in rows if r[0].startswith("data/delivered") and r[4]))
