"""Release validation: CURIEs against source term sets, gene symbol/id consistency, Biolink categories and predicates
(bmt, Biolink 4.4.6), dangling edges, study-result links, and the KGX validator (streaming)."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from .biolink import toolkit, visionkb_classes
from .ids import check_curie, prefix
from .kgx import read_jsonl

KNOWLEDGE_LEVELS = {
    "knowledge_assertion",
    "logical_entailment",
    "prediction",
    "statistical_association",
    "text_co_occurrence",
    "observation",
    "not_provided",
}
AGENT_TYPES = {
    "manual_agent",
    "automated_agent",
    "data_analysis_pipeline",
    "computational_model",
    "text_mining_agent",
    "image_processing_agent",
    "manual_validation_of_automated_agent",
    "not_provided",
}
MAX_EXAMPLES = 20


class Report:
    def __init__(self):
        self.counts = Counter()
        self.examples = defaultdict(list)
        self.level = {}

    def add(self, level: str, code: str, example):
        self.counts[code] += 1
        self.level[code] = level
        if len(self.examples[code]) < MAX_EXAMPLES:
            self.examples[code].append(example)

    def summary(self) -> dict:
        errs = {k: v for k, v in self.counts.items() if self.level[k] == "error"}
        warns = {k: v for k, v in self.counts.items() if self.level[k] == "warning"}
        return {
            "n_errors": sum(errs.values()),
            "n_warnings": sum(warns.values()),
            "errors": {k: {"count": v, "examples": self.examples[k]} for k, v in sorted(errs.items())},
            "warnings": {k: {"count": v, "examples": self.examples[k]} for k, v in sorted(warns.items())},
        }


class TermSets:
    """Valid ids per prefix. hgnc: id -> approved symbol; obo: prefix -> set of non-obsolete ids."""

    def __init__(self, hgnc: dict[str, str], obo: dict[str, set[str]]):
        self.hgnc = hgnc
        self.obo = obo


def check_node(n: dict, ts: TermSets, rep: Report, vk_classes: dict, tk):
    nid = n.get("id")
    why = check_curie(nid)
    if why:
        rep.add("error", "node_id_malformed", {"id": nid, "reason": why})
        return
    cats = n.get("category") or []
    if not cats:
        rep.add("error", "node_without_category", nid)
    for c in cats:
        if c.startswith("biolink:"):
            if not tk.is_category(c):
                rep.add("error", "unknown_biolink_category", {"id": nid, "category": c})
        elif c.startswith("visionkb:"):
            if c.split(":", 1)[1] not in vk_classes:
                rep.add("error", "unknown_visionkb_category", {"id": nid, "category": c})
        else:
            rep.add("error", "category_not_curie", {"id": nid, "category": c})
    p = prefix(nid)
    if p == "HGNC":
        sym = ts.hgnc.get(nid)
        if sym is None:
            rep.add("error", "hgnc_id_not_approved", nid)
        elif n.get("symbol") and n["symbol"] != sym:
            rep.add("error", "hgnc_symbol_mismatch", {"id": nid, "symbol": n["symbol"], "hgnc_symbol": sym})
    elif p in ts.obo and nid not in ts.obo[p]:
        rep.add("error", "ontology_id_unknown_or_obsolete", nid)
    elif p in ("EFO", "NCBITaxon") and not n.get("ols_verified"):
        rep.add("error", "ontology_id_not_verified", nid)
    if not n.get("name"):
        rep.add("warning", "node_without_name", nid)
    if "biolink:StudyResult" in cats:
        t = n.get("supporting_text") or []
        for x in [t] if isinstance(t, str) else t:
            if len(x) > 300:
                rep.add("error", "supporting_text_too_long", {"id": nid, "len": len(x)})


def validate_release(nodes_path: Path, edges_path: Path, ts: TermSets, *, run_kgx: bool = True, log=print) -> dict:
    tk = toolkit()
    vk = visionkb_classes()
    rep = Report()
    node_cat: dict[str, tuple] = {}
    sr_edge: dict[str, str] = {}
    datasets = set()
    by_cat = Counter()
    log("validating nodes")
    for n in read_jsonl(nodes_path):
        check_node(n, ts, rep, vk, tk)
        cats = tuple(n.get("category") or ())
        node_cat[n["id"]] = cats
        for c in cats:
            by_cat[c] += 1
        if "biolink:StudyResult" in cats:
            sr_edge[n["id"]] = n.get("edge_id")
            for d in n.get("supporting_data_set") or []:
                datasets.add((n["id"], d))
        if "biolink:Dataset" in cats:
            pass
    for sid, d in datasets:
        if d not in node_cat:
            rep.add("error", "study_result_dataset_missing", {"study_result": sid, "dataset": d})
    log("validating edges")
    by_pred = Counter()
    triples = Counter()
    edge_ids = set()
    sr_used = Counter()
    for e in read_jsonl(edges_path):
        eid, s, p, o = e.get("id"), e.get("subject"), e.get("predicate"), e.get("object")
        by_pred[p] += 1
        if eid in edge_ids:
            rep.add("error", "duplicate_edge_id", eid)
        edge_ids.add(eid)
        for end, v in (("subject", s), ("object", o)):
            if v not in node_cat:
                rep.add("error", f"dangling_{end}", {"edge": eid, end: v})
        if not tk.is_predicate(p):
            rep.add("error", "unknown_predicate", {"edge": eid, "predicate": p})
        if e.get("knowledge_level") not in KNOWLEDGE_LEVELS:
            rep.add("error", "bad_knowledge_level", {"edge": eid, "value": e.get("knowledge_level")})
        if e.get("agent_type") not in AGENT_TYPES:
            rep.add("error", "bad_agent_type", {"edge": eid, "value": e.get("agent_type")})
        if check_curie(e.get("primary_knowledge_source") or "") or not str(
            e.get("primary_knowledge_source")
        ).startswith("infores:"):
            rep.add("error", "bad_primary_knowledge_source", eid)
        for q in e.get("anatomical_context_qualifier") or []:
            if q not in node_cat:
                rep.add("error", "qualifier_term_not_in_graph", {"edge": eid, "value": q})
        sq = e.get("species_context_qualifier")
        if sq and sq not in node_cat:
            rep.add("error", "qualifier_term_not_in_graph", {"edge": eid, "value": sq})
        for pm in e.get("publications") or []:
            if check_curie(pm) or not pm.startswith("PMID:"):
                rep.add("error", "bad_publication_id", {"edge": eid, "value": pm})
        for sid in e.get("has_study_results") or []:
            sr_used[sid] += 1
            if sid not in sr_edge:
                rep.add("error", "study_result_missing", {"edge": eid, "study_result": sid})
            elif sr_edge[sid] != eid:
                rep.add("error", "study_result_edge_id_mismatch", {"edge": eid, "study_result": sid})
        if s in node_cat and o in node_cat:
            triples[(node_cat[s][0], p, node_cat[o][0])] += 1
    for sid in sr_edge:
        if sid not in sr_used:
            rep.add("error", "study_result_orphan", sid)
    domain_range = {}
    for (sc, p, oc), n in sorted(triples.items()):
        sc_b = (
            sc
            if sc.startswith("biolink:")
            else "biolink:" + vk.get(sc.split(":")[1], "named thing").title().replace(" ", "")
        )
        oc_b = (
            oc
            if oc.startswith("biolink:")
            else "biolink:" + vk.get(oc.split(":")[1], "named thing").title().replace(" ", "")
        )
        ok = tk.validate_edge(sc_b, p, oc_b)
        domain_range[f"{sc} {p} {oc}"] = {"count": n, "bmt_valid": bool(ok)}
        if not ok:
            rep.add("warning", "bmt_domain_range_not_valid", {"triple": f"{sc} {p} {oc}", "count": n})
    out = {
        "biolink_version": tk.get_model_version(),
        "nodes_by_category": dict(sorted(by_cat.items())),
        "edges_by_predicate": dict(sorted(by_pred.items())),
        "category_predicate_triples": domain_range,
    }
    if run_kgx:
        log("running KGX validator (streaming)")
        out["kgx"] = kgx_validate(nodes_path, edges_path, tk)
        for code, info in out["kgx"]["expected"].items():
            rep.counts[f"kgx:expected:{code}"] += info["count"]
            rep.level[f"kgx:expected:{code}"] = "warning"
            rep.examples[f"kgx:expected:{code}"] = info["examples"]
        for lvl, d in out["kgx"]["errors"].items():
            for code, info in d.items():
                rep.counts[f"kgx:{lvl}:{code}"] += info["count"]
                rep.level[f"kgx:{lvl}:{code}"] = "error" if lvl == "ERROR" else "warning"
                rep.examples[f"kgx:{lvl}:{code}"] = info["examples"]
    out.update(rep.summary())
    out["status"] = "pass" if out["n_errors"] == 0 else "fail"
    return out


def kgx_validate(nodes_path: Path, edges_path: Path, tk) -> dict:
    """KGX 2.8 validator (streaming) with the vendored Biolink 4.4.6 toolkit.

    Expected findings are reported as warnings, not hidden: the local ``visionkb`` prefix (CellState, Dataset, StudyResult ids;
    defined in the visionkb schema, not in the Biolink JSON-LD context), visionkb categories (``CellState``), and a KGX streaming
    artefact where the last node before the edge stream is re-checked with list-valued ``in_taxon`` (each such node is re-read
    from the file and downgraded only if its ``in_taxon`` is a string there). Everything else stays an error.
    """
    from kgx.transformer import Transformer
    from kgx.validator import Validator

    Validator._currently_active_toolkit = tk
    v = Validator()
    t = Transformer(stream=True)
    t.transform(
        input_args={"filename": [str(nodes_path), str(edges_path)], "format": "jsonl"},
        output_args={"format": "null"},
        inspector=v,
    )
    errs = v.get_errors()
    vk = {f"'{c}'" for c in visionkb_classes()}
    taxon_ids = set()
    for codes in errs.values():
        for msgs in codes.values():
            for m, ids in msgs.items():
                if "in_taxon" in m and "Single-valued" in m:
                    taxon_ids |= set(ids if isinstance(ids, list | set) else [ids])
    taxon_ok = set()
    if taxon_ids:
        for n in read_jsonl(nodes_path):
            if n["id"] in taxon_ids and isinstance(n.get("in_taxon"), str):
                taxon_ok.add(n["id"])
    summary = {}
    expected = {}
    for lvl, codes in errs.items():
        for code, msgs in codes.items():
            for m, ids in msgs.items():
                idl = sorted(set(ids)) if isinstance(ids, list | set) else [ids]
                if "CURIE prefix 'visionkb'" in m:
                    bucket, key = expected, f"{code}:local_visionkb_prefix"
                elif code == "INVALID_CATEGORY" and any(c in m for c in vk):
                    bucket, key = expected, f"{code}:visionkb_category"
                elif "in_taxon" in m and idl and set(idl) <= taxon_ok:
                    bucket, key = expected, f"{code}:kgx_stream_artifact_in_taxon"
                else:
                    bucket, key = summary.setdefault(lvl, {}), code
                rec = bucket.setdefault(key, {"count": 0, "examples": []})
                rec["count"] += len(idl)
                if len(rec["examples"]) < 5:
                    rec["examples"].append({"message": m[:300], "ids": idl[:5]})
    return {"errors": summary, "expected": expected, "biolink_version": tk.get_model_version()}


def apln_self_test(ts: TermSets) -> dict:
    """The design proposal used HGNC:16693 for APLN (APLN is HGNC:16665; HGNC:16693 is another gene): must be flagged."""
    tk = toolkit()
    vk = visionkb_classes()
    rep = Report()
    probe = {"id": "HGNC:16693", "category": ["biolink:Gene"], "name": "APLN", "symbol": "APLN"}
    check_node(probe, ts, rep, vk, tk)
    good = Report()
    check_node({"id": "HGNC:16665", "category": ["biolink:Gene"], "name": "APLN", "symbol": "APLN"}, ts, good, vk, tk)
    flagged = rep.counts.get("hgnc_symbol_mismatch", 0) > 0
    return {
        "probe": probe,
        "flagged": flagged,
        "found": rep.summary()["errors"],
        "hgnc_symbol_of_HGNC:16693": ts.hgnc.get("HGNC:16693"),
        "correct_id_passes": good.summary()["n_errors"] == 0,
        "status": "pass" if flagged and good.summary()["n_errors"] == 0 else "fail",
    }
