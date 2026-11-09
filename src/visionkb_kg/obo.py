"""Minimal OBO 1.4 reader: terms with names, synonyms, xrefs, is_a and part_of, obsolescence."""

from __future__ import annotations

import re
from pathlib import Path

SYN = re.compile(r'^"((?:[^"\\]|\\.)*)"\s+(EXACT|BROAD|NARROW|RELATED)')


def read_obo(path: str | Path, prefixes: tuple[str, ...] | None = None) -> tuple[dict, dict]:
    """Return (header, terms). terms[id] = {name, def, synonyms, exact_synonyms, xrefs, is_a, part_of, obsolete, replaced_by}."""
    header, terms = {}, {}
    cur = None
    in_header = True
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_header = False
                cur = {} if line == "[Term]" else None
                if cur is not None:
                    cur.update(
                        synonyms=[], exact_synonyms=[], xrefs=[], is_a=[], part_of=[], obsolete=False, replaced_by=[]
                    )
                continue
            if not line or line.startswith("!"):
                if not line and cur is not None and cur.get("id"):
                    terms[cur["id"]] = cur
                    cur = None
                continue
            k, _, v = line.partition(": ")
            if in_header:
                header.setdefault(k, v)
                continue
            if cur is None:
                continue
            v = re.sub(r"\s+!.*$", "", v) if k in ("is_a", "relationship", "replaced_by") else v
            if k == "id":
                cur["id"] = v
            elif k == "name":
                cur["name"] = v
            elif k == "def":
                m = re.match(r'^"((?:[^"\\]|\\.)*)"', v)
                cur["def"] = m.group(1).replace('\\"', '"') if m else v
            elif k == "synonym":
                m = SYN.match(v)
                if m:
                    s = m.group(1).replace('\\"', '"')
                    cur["synonyms"].append(s)
                    if m.group(2) == "EXACT":
                        cur["exact_synonyms"].append(s)
            elif k == "xref":
                cur["xrefs"].append(v.split(" ")[0])
            elif k == "is_a":
                cur["is_a"].append(v.split(" ")[0])
            elif k == "relationship":
                parts = v.split()
                if len(parts) >= 2 and parts[0] in ("part_of", "BFO:0000050"):
                    cur["part_of"].append(parts[1])
            elif k == "is_obsolete":
                cur["obsolete"] = v.strip() == "true"
            elif k == "replaced_by":
                cur["replaced_by"].append(v.strip())
    if cur is not None and cur.get("id"):
        terms[cur["id"]] = cur
    if prefixes:
        terms = {k: v for k, v in terms.items() if k.split(":")[0] in prefixes}
    return header, terms


def version_of(header: dict) -> str:
    return header.get("data-version") or header.get("date") or "unknown"


def ancestors(terms: dict, seeds, rels=("is_a", "part_of")) -> set:
    out, todo = set(), list(seeds)
    while todo:
        t = todo.pop()
        if t in out:
            continue
        out.add(t)
        rec = terms.get(t)
        if rec:
            for r in rels:
                todo.extend(rec[r])
    return out


def descendants(terms: dict, root: str, rels=("is_a", "part_of")) -> set:
    children: dict[str, list[str]] = {}
    for tid, rec in terms.items():
        if rec["obsolete"]:
            continue
        for r in rels:
            for p in rec[r]:
                children.setdefault(p, []).append(tid)
    out, todo = set(), [root]
    while todo:
        t = todo.pop()
        if t in out:
            continue
        out.add(t)
        todo.extend(children.get(t, []))
    return out
