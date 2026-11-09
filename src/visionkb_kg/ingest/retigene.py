"""Retigene curated haploinsufficiency review -> gene to disease (MONDO) / phenotype (HP) edges with publications.

Grades (curator review_grade): A -> biolink:causes; B (hedged/conflicting) -> biolink:gene_associated_with_condition
(MONDO) or biolink:has_phenotype (HP);
N -> no edge (reported). subject_form_or_variant_qualifier = loss_of_function_variant_form (heterozygous LoF). Phenotype
text is mapped through a curated query table searched in OLS: a term is accepted only when its label or an exact
synonym equals the query; nothing is guessed, unmapped genes are reported. Evidence sentences become
LiteratureStudyResult nodes (supporting_text: Biolink multivalued, one sentence <= 300 characters).
"""

from __future__ import annotations

import time

import pandas as pd
import requests

from ..ids import short_id
from .common import make_edge

SOURCE = "infores:visionkb-retigene"
OLS = "https://www.ebi.ac.uk/ols4/api/search"
# A -> causes; B (hedged/conflicting) -> gene_associated_with_condition (MONDO) or has_phenotype (HP)
PREDICATE = {"A": "biolink:causes", "B": "biolink:gene_associated_with_condition"}


def predicate_for(grade: str, obj: str) -> str:
    if grade == "B" and obj.startswith("HP:"):
        return "biolink:has_phenotype"
    return PREDICATE[grade]


MAX_TEXT = 300


def clip(text, n=MAX_TEXT) -> str | None:
    if not isinstance(text, str) or not text.strip():
        return None
    t = " ".join(text.split())
    return t if len(t) <= n else t[: n - 1].rstrip() + "\u2026"


def ols_exact(query: str, ontology: str, session=None) -> list[dict]:
    """Terms of `ontology` whose label or exact synonym equals query (case-insensitive), not obsolete."""
    s = session or requests
    pre = {"mondo": "MONDO", "hp": "HP"}[ontology]
    for attempt in range(5):
        try:
            r = s.get(
                OLS,
                params={
                    "q": query,
                    "ontology": ontology,
                    "rows": 25,
                    "queryFields": "label,synonym",
                    "exact": "true",
                    "fieldList": "obo_id,label,synonym,exact_synonyms,is_obsolete",
                },
                timeout=60,
            )
            r.raise_for_status()
            docs = r.json()["response"]["docs"]
            break
        except (requests.RequestException, KeyError, ValueError):
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))
    q = query.casefold()
    out = []
    for d in docs:
        oid = d.get("obo_id") or ""
        if not oid.startswith(pre + ":") or d.get("is_obsolete"):
            continue
        syns = [x.casefold() for x in (d.get("exact_synonyms") or d.get("synonym") or [])]
        if d.get("label", "").casefold() == q:
            out.append({"id": oid, "label": d["label"], "match": "label"})
        elif q in syns:
            out.append(
                {"id": oid, "label": d["label"], "match": "exact_synonym" if d.get("exact_synonyms") else "synonym"}
            )
    return out


def map_phenotypes(genes: pd.DataFrame, queries: dict) -> pd.DataFrame:
    rows = []
    sess = requests.Session()
    for r in genes.itertuples(index=False):
        qs = queries.get(r.gene)
        if qs is None:
            rows.append({"gene": r.gene, "phenotype": r.phenotype, "status": "no_query", "object_id": None})
            continue
        if not qs:
            rows.append(
                {"gene": r.gene, "phenotype": r.phenotype, "status": "unmapped_no_named_condition", "object_id": None}
            )
            continue
        got = False
        for q in qs:
            hits = ols_exact(q["query"], q["ontology"], sess)
            if len(hits) == 1:
                h = hits[0]
                rows.append(
                    {
                        "gene": r.gene,
                        "phenotype": r.phenotype,
                        "query": q["query"],
                        "status": "mapped",
                        "object_id": h["id"],
                        "object_label": h["label"],
                        "match": h["match"],
                        "confidence": q["confidence"],
                        "note": q.get("note"),
                    }
                )
                got = True
            elif len(hits) > 1:
                rows.append(
                    {
                        "gene": r.gene,
                        "phenotype": r.phenotype,
                        "query": q["query"],
                        "status": "ambiguous",
                        "object_id": None,
                        "candidates": ";".join(h["id"] for h in hits),
                        "confidence": q["confidence"],
                    }
                )
            else:
                rows.append(
                    {
                        "gene": r.gene,
                        "phenotype": r.phenotype,
                        "query": q["query"],
                        "status": "no_exact_match",
                        "object_id": None,
                        "confidence": q["confidence"],
                    }
                )
        if not got:
            pass
    df = pd.DataFrame(rows)
    # one row per (gene, object): first accepted query wins
    return df


def build(genes: pd.DataFrame, evidence: pd.DataFrame, mapping: pd.DataFrame, hgnc_symbols: dict, release: str):
    nodes, edges, flagged = [], [], []
    for r in genes.sort_values("gene").itertuples(index=False):
        grade = str(r.review_grade).strip()
        hid = hgnc_symbols.get(r.gene)
        if hid is None:
            flagged.append({"gene": r.gene, "reason": "symbol not an approved HGNC symbol"})
            continue
        if grade not in PREDICATE:
            flagged.append({"gene": r.gene, "grade": grade, "reason": "grade N (not supported): no edge"})
            continue
        objs = mapping[(mapping.gene == r.gene) & (mapping.status == "mapped")].drop_duplicates("object_id")
        if objs.empty:
            st = ";".join(sorted(set(mapping.loc[mapping.gene == r.gene, "status"]))) or "no_query"
            flagged.append(
                {
                    "gene": r.gene,
                    "grade": grade,
                    "phenotype": r.phenotype,
                    "reason": f"phenotype unmapped ({st}): no edge",
                }
            )
            continue
        ev = evidence[evidence.gene == r.gene]
        dosage = ev[ev["dosage_wording"].fillna("").astype(str).str.strip() != ""]
        pmids = sorted(
            {f"PMID:{int(p)}" for p in [r.pmid, *dosage["pmid"].tolist()] if pd.notna(p)}, key=lambda x: int(x[5:])
        )
        for o in objs.itertuples(index=False):
            e = make_edge(
                hid,
                predicate_for(grade, o.object_id),
                o.object_id,
                source=SOURCE,
                knowledge_level="knowledge_assertion",
                agent_type="manual_validation_of_automated_agent",
                release=release,
                subject_form_or_variant_qualifier="loss_of_function_variant_form",
                publications=pmids,
                grade=grade,
                status="curated",
                visibility="lab",
                phenotype_text=r.phenotype,
                mapping_confidence=o.confidence,
                mapping_match=o.match,
            )
            srs = []
            # cited reference sentence first, then sentences with dosage wording
            ref = ev[ev["is_reference_quote"].astype(str) == "yes"]
            seen = set()
            for kind, sub in (("reference_quote", ref), ("dosage_sentence", dosage)):
                for s in sub.itertuples(index=False):
                    key = (s.pmid, s.sentence_id)
                    if key in seen:
                        continue
                    seen.add(key)
                    sid = short_id("sr", e["id"], s.pmid, s.sentence_id)
                    srs.append(sid)
                    nodes.append(
                        {
                            "id": sid,
                            "category": ["biolink:StudyResult"],
                            "name": f"{r.gene} haploinsufficiency evidence (PMID:{int(s.pmid)})",
                            "edge_id": e["id"],
                            "pmid": f"PMID:{int(s.pmid)}",
                            "supporting_text": [clip(s.quote)] if clip(s.quote) else None,
                            "grade": grade,
                            "review_note": clip(r.review_note, 1000) if kind == "reference_quote" else None,
                            "evidence_kind": kind,
                            "page": int(s.page) if pd.notna(s.page) else None,
                            "provided_by": [SOURCE],
                        }
                    )
            e["has_study_results"] = sorted(srs)
            e["evidence_count"] = len(srs)
            edges.append(e)
    nodes = list({n["id"]: n for n in nodes}.values())
    stats = {
        "genes": len(genes),
        "edges": len(edges),
        "study_results": len(nodes),
        "edges_by_predicate": pd.Series([e["predicate"] for e in edges]).value_counts().to_dict() if edges else {},
        "flagged": flagged,
    }
    return nodes, edges, stats
