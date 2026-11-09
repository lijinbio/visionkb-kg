"""Build mappings/heca-celltypes.sssom.tsv: atlas label -> Cell Ontology, OLS-checked.

Input: a JSON file {dataset_id: {"majorclass": [labels], "celltype": {label: dominant majorclass}}} listing the labels
present in each atlas (written from the atlas obs; kept outside git). Output: one SSSOM row per (dataset, level, label)
that has a curated mapping; labels without a row become CellStates (contracts §3).

Rule (reused from the 2026-11-05 prototype): the CL term encodes cell-type identity only; tissue context is carried by
the edge's anatomical_context_qualifier. exactMatch = the label names that CL class; broadMatch = the CL class is
broader than the label but is the best class for it; closeMatch = similar, no adequate CL class (the label becomes a
CellState subclass of it). Every CL id is checked against OLS (label equal to the expected label, not obsolete).

Usage: python scripts/build_crosswalk.py labels.json [--no-ols]
"""

from __future__ import annotations

import csv
import datetime
import json
import re
import sys
import time
from pathlib import Path

import requests

OUT = Path(__file__).resolve().parents[1] / "src" / "visionkb_kg" / "mappings" / "heca-celltypes.sssom.tsv"

E, B, C = "skos:exactMatch", "skos:broadMatch", "skos:closeMatch"
# majorclass labels (prototype s02 crosswalk, reviewed 2026-11-05)
M = {
    "AC": ("CL:0000561", "amacrine cell", E, 0.95, ""),
    "BC": ("CL:0000748", "retinal bipolar neuron", E, 0.95, ""),
    "Cone": ("CL:0000573", "retinal cone cell", E, 0.95, ""),
    "HC": ("CL:0000745", "retina horizontal cell", E, 0.95, ""),
    "MG": ("CL:0000636", "Mueller cell", E, 0.95, ""),
    "Microglia": ("CL:0000129", "microglial cell", E, 0.95, ""),
    "RGC": ("CL:0000740", "retinal ganglion cell", E, 0.95, ""),
    "Rod": ("CL:0000604", "retinal rod cell", E, 0.95, ""),
    "RPE": ("CL:0002586", "retinal pigment epithelial cell", E, 0.95, ""),
    "RPE cells": ("CL:0002586", "retinal pigment epithelial cell", E, 0.95, ""),
    "RPE_Fetal": (
        "CL:0002586",
        "retinal pigment epithelial cell",
        E,
        0.9,
        "fetal stage belongs on the donor, not in the cell-type term",
    ),
    "Anterior Pigmented Epithelial": (
        "CL:0002565",
        "iris pigment epithelial cell",
        B,
        0.6,
        "anterior layer is the dilator myoepithelium; candidate CL:1000445",
    ),
    "Posterior Pigmented Epithelial": ("CL:0002565", "iris pigment epithelial cell", C, 0.8, ""),
    "Putative Stem Cell": ("CL:0000034", "stem cell", B, 0.4, "putative population; needs expert review"),
    "Sphincter Muscle": ("CL:0002243", "smooth muscle cell of sphincter of pupil", E, 0.9, ""),
    "Smooth Muscle": ("CL:0000192", "smooth muscle cell", E, 0.9, ""),
    "Smooth Muscle Cells": ("CL:0000192", "smooth muscle cell", E, 0.9, ""),
    "Vascular Endothelial": ("CL:0000071", "blood vessel endothelial cell", E, 0.9, ""),
    "Vascular_endothelium": ("CL:0000071", "blood vessel endothelial cell", E, 0.9, ""),
    "Anterior_Epi": ("CL:0002223", "anterior lens cell", E, 0.85, ""),
    "Equatorial_Epi": ("CL:0002224", "lens epithelial cell", B, 0.8, "equatorial (germinative) lens epithelium"),
    "Early_Fiber": ("CL:0011004", "lens fiber cell", B, 0.7, "newly differentiating fibers"),
    "Fiber": ("CL:0011004", "lens fiber cell", E, 0.9, ""),
    "Transition": (
        "CL:0002224",
        "lens epithelial cell",
        C,
        0.4,
        "epithelial-to-fiber transition zone; no CL class; term-request candidate",
    ),
    "Corneal Endothelium": ("CL:0000132", "corneal endothelial cell", E, 0.95, ""),
    "Epithelium": (
        "CL:0000066",
        "epithelial cell",
        E,
        0.7,
        "heterogeneous (corneal, limbal, conjunctival); resolve at subclass level",
    ),
    "Goblet": ("CL:2000084", "conjunctiva goblet cell", E, 0.9, ""),
    "Ciliary_Muscle": ("CL:1000443", "ciliary muscle cell", E, 0.9, ""),
    "CB_NPCE": ("CL:0002304", "non-pigmented ciliary epithelial cell", E, 0.95, ""),
    "CB_PCE": ("CL:0002303", "pigmented ciliary epithelial cell", E, 0.95, ""),
    "Oligodendrocyte": ("CL:0000128", "oligodendrocyte", E, 0.95, ""),
    "Oligodendrocyte_precursor_cell": ("CL:0002453", "oligodendrocyte precursor cell", E, 0.95, ""),
    "OPC": ("CL:0002453", "oligodendrocyte precursor cell", E, 0.95, ""),
    "Lymphocyte": ("CL:0000542", "lymphocyte", E, 0.9, ""),
    "DC": ("CL:0000451", "dendritic cell", E, 0.9, ""),
    "pDC": ("CL:0000784", "plasmacytoid dendritic cell", E, 0.95, ""),
    "ILC3": ("CL:0001071", "group 3 innate lymphoid cell", E, 0.95, ""),
    "Megakaryocytes": ("CL:0000556", "megakaryocyte", E, 0.9, ""),
    "Astrocyte": ("CL:0000127", "astrocyte", E, 0.9, ""),
}
_GROUPS = {
    ("CL:0000236", "B cell", 0.95): ["B_cells", "B cells", "B Cells", "B_cell", "B Cell"],
    ("CL:0000084", "T cell", 0.95): ["T_cells", "T cells", "T_cell"],
    ("CL:0000623", "natural killer cell", 0.95): ["NK_cells", "NK cells", "NK_cell"],
    ("CL:0000786", "plasma cell", 0.95): ["Plasma_cells", "Plasma cells", "Plasma Cells"],
    ("CL:0000097", "mast cell", 0.95): ["Mast_cells", "Mast cells", "Mast Cells", "Mast_cell", "Mast Cell"],
    ("CL:0000451", "dendritic cell", 0.95): ["Dendritic_cells", "Dendritic_cell"],
    ("CL:0000057", "fibroblast", 0.9): ["Fibroblast", "Fibroblasts"],
    ("CL:0000148", "melanocyte", 0.9): ["Melanocyte", "Melanocytes"],
    ("CL:0000576", "monocyte", 0.95): ["Monocyte", "Monocytes"],
    ("CL:0000235", "macrophage", 0.95): ["Macrophage", "Macrophages"],
    ("CL:0002573", "Schwann cell", 0.95): [
        "Schwann",
        "Schwann Cells",
        "Schwann_Cell",
        "Schwann Cell",
        "Schwann_cells",
        "Schwann cells",
    ],
    ("CL:0008034", "mural cell", 0.9): ["Mural Cells", "Mural_cell", "Mural cells", "Mural_cells", "Mural"],
    ("CL:0000669", "pericyte", 0.95): ["Pericyte", "Pericytes"],
    ("CL:0000115", "endothelial cell", 0.85): ["Endothelium", "Endothelial_cell", "Endothelial cells", "Endothelial"],
    ("CL:0000738", "leukocyte", 0.8): ["Immune_Cell", "Immune Cell", "Immune Cells", "Immune"],
}
for (cid, lab, conf), labels in _GROUPS.items():
    for x in labels:
        M.setdefault(x, (cid, lab, E, conf, "mixed immune population" if lab == "leukocyte" else ""))

# cell-type level labels with a CL class (checked in OLS 2026-11-09); other subtypes become CellStates under their majorclass
CT = {
    "Rod": M["Rod"],
    "RPE": M["RPE"],
    "MG": M["MG"],
    "Microglia": M["Microglia"],
    "Goblet": M["Goblet"],
    "OPC": M["OPC"],
    "DC": M["DC"],
    "pDC": M["pDC"],
    "Megakaryocytes": M["Megakaryocytes"],
    "CB_PCE": M["CB_PCE"],
    "Ciliary_Muscle": M["Ciliary_Muscle"],
    "Smooth Muscle Cells": M["Smooth Muscle Cells"],
    "S_Cone": ("CL:0003050", "S cone cell", E, 0.95, ""),
    "RB": ("CL:0000751", "rod bipolar cell", E, 0.9, ""),
    "DB1": ("CL:4033027", "diffuse bipolar 1 cell", E, 0.9, ""),
    "DB2": ("CL:4033028", "diffuse bipolar 2 cell", E, 0.9, ""),
    "DB3a": ("CL:4033029", "diffuse bipolar 3a cell", E, 0.9, ""),
    "DB3b": ("CL:4033030", "diffuse bipolar 3b cell", E, 0.9, ""),
    "DB4a": ("CL:4033095", "diffuse bipolar 4a cell", E, 0.9, ""),
    "DB4b": ("CL:4033096", "diffuse bipolar 4b cell", E, 0.9, ""),
    "DB5": ("CL:4033085", "diffuse bipolar 5 cell", E, 0.9, ""),
    "DB6": ("CL:4033032", "diffuse bipolar 6 cell", E, 0.9, ""),
    "FMB": ("CL:4033033", "flat midget bipolar cell", E, 0.9, ""),
    "IMB": ("CL:4033034", "invaginating midget bipolar cell", E, 0.9, ""),
    "GB": ("CL:4033035", "giant bipolar cell", E, 0.85, ""),
    "BB": ("CL:4033019", "ON-blue cone bipolar cell", E, 0.85, "BB = blue-cone bipolar"),
    "H1": ("CL:0004217", "H1 horizontal cell", E, 0.9, ""),
    "H2": ("CL:0004218", "H2 horizontal cell", E, 0.9, ""),
    "ipRGC": ("CL:0020014", "intrinsically photosensitive retinal ganglion cell", E, 0.9, ""),
    "MG_ON": ("CL:4033046", "ON midget ganglion cell", E, 0.9, "MG here is midget ganglion (RGC), not Mueller glia"),
    "MG_OFF": ("CL:4033047", "OFF midget ganglion cell", E, 0.9, "MG here is midget ganglion (RGC), not Mueller glia"),
    "PG_ON": ("CL:4033052", "ON parasol ganglion cell", E, 0.9, ""),
    "PG_OFF": ("CL:4033051", "OFF parasol ganglion cell", E, 0.9, ""),
    "Neutrophils": ("CL:0000775", "neutrophil", E, 0.9, ""),
    "Corneal_Endo": M["Corneal Endothelium"],
    "Lymphatic_Endo": ("CL:0002138", "endothelial cell of lymphatic vessel", E, 0.9, ""),
    "Lymphatic Endothelium": ("CL:0002138", "endothelial cell of lymphatic vessel", E, 0.9, ""),
    "VSMC": ("CL:0000359", "vascular associated smooth muscle cell", E, 0.9, ""),
    "Vascular smooth muscle cells": ("CL:0000359", "vascular associated smooth muscle cell", E, 0.9, ""),
    "Smooth_muscle_cells": M["Smooth Muscle"],
    "Schlemm Endothelium": ("CL:4033097", "Schlemm's canal endothelial cell", E, 0.9, ""),
    "Schlemm_Endothelium": ("CL:4033097", "Schlemm's canal endothelial cell", E, 0.9, ""),
    "Vascular Endothelium": M["Vascular Endothelial"],
    "Vascular_Endothelium": M["Vascular Endothelial"],
    "Artery": ("CL:1000413", "endothelial cell of artery", E, 0.8, "endothelial subcluster labelled by vessel"),
    "Arterial": ("CL:1000413", "endothelial cell of artery", E, 0.8, "endothelial subcluster labelled by vessel"),
    "Capillary": ("CL:0002144", "capillary endothelial cell", E, 0.85, "endothelial subcluster labelled by vessel"),
    "Capillary_fenestrated": ("CL:0000666", "fenestrated endothelial cell", E, 0.8, ""),
    "Fenestrated_capillary": ("CL:0000666", "fenestrated endothelial cell", E, 0.8, ""),
    "Fenestrated capillary": ("CL:0000666", "fenestrated endothelial cell", E, 0.8, ""),
    "Vein": ("CL:0002543", "vein endothelial cell", E, 0.85, ""),
    "Venule": ("CL:1000414", "endothelial cell of venule", E, 0.85, ""),
    "Venules": ("CL:1000414", "endothelial cell of venule", E, 0.85, ""),
    "Corneal Stromal Keratocytes": ("CL:0002363", "keratocyte", E, 0.9, ""),
    "Conj-Epi": ("CL:1000432", "conjunctival epithelial cell", E, 0.85, ""),
    "T Cells_CD4": ("CL:0000624", "CD4-positive, alpha-beta T cell", E, 0.8, ""),
    "T Cells_CD8": ("CL:0000625", "CD8-positive, alpha-beta T cell", E, 0.8, ""),
    "Helper T cells": ("CL:0000912", "helper T cell", E, 0.9, ""),
    "Regulatory T cells": ("CL:0000815", "regulatory T cell", E, 0.9, ""),
    "Classic Monocytes": ("CL:0000860", "classical monocyte", E, 0.9, ""),
    "Non-classical Monocytes": ("CL:0000875", "non-classical monocyte", E, 0.9, ""),
    "ILC": ("CL:0001065", "innate lymphoid cell", E, 0.9, ""),
    "Schwann_M": ("CL:0000218", "myelinating Schwann cell", E, 0.9, ""),
    "Schwann-my": ("CL:0000218", "myelinating Schwann cell", E, 0.9, ""),
    "Schwann_myelinating": ("CL:0000218", "myelinating Schwann cell", E, 0.9, ""),
    "Schwann cells myelinating": ("CL:0000218", "myelinating Schwann cell", E, 0.9, ""),
    "Schwann_N": ("CL:0002376", "non-myelinating Schwann cell", E, 0.9, ""),
    "Schwann-nmy": ("CL:0002376", "non-myelinating Schwann cell", E, 0.9, ""),
    "Schwann_nonmyelinating": ("CL:0002376", "non-myelinating Schwann cell", E, 0.9, ""),
    "Schwann cells nonmyelinating": ("CL:0002376", "non-myelinating Schwann cell", E, 0.9, ""),
    "TM_Fibro_JCT": ("CL:7770002", "juxtacanalicular tissue cell", B, 0.6, "JCT fibroblast cluster"),
    "JCT": ("CL:7770002", "juxtacanalicular tissue cell", B, 0.6, "JCT fibroblast cluster"),
    "Astro_retina": ("CL:4033015", "retinal astrocyte", E, 0.85, ""),
}
for lab in [
    "B_cells",
    "B cells",
    "B Cells",
    "B_cell",
    "B Cell",
    "T_cells",
    "T_cell",
    "NK_cell",
    "Plasma_cells",
    "Plasma cells",
    "Plasma Cells",
    "Mast_cells",
    "Mast cells",
    "Mast Cells",
    "Mast_cell",
    "Mast Cell",
    "Dendritic_cells",
    "Dendritic_cell",
    "Melanocyte",
    "Melanocytes",
    "Monocytes",
    "Macrophage",
    "Macrophages",
    "Pericytes",
    "Astrocyte",
]:
    CT.setdefault(lab, M[lab])
# (dataset_id, label) overrides at any level
OVR = {
    ("visionkb:ds-hrca-scrna", "Astrocyte"): ("CL:4033015", "retinal astrocyte", E, 0.9, ""),
    ("visionkb:ds-hrca-snrna", "Astrocyte"): ("CL:4033015", "retinal astrocyte", E, 0.9, ""),
    ("visionkb:ds-ononh-snrna", "Astrocyte"): (
        "CL:0000127",
        "astrocyte",
        E,
        0.9,
        "optic-nerve astrocytes; no tissue-specific CL class",
    ),
    ("visionkb:ds-rpechoroid-scrna", "Melanocytes"): ("CL:4030000", "choroidal melanocyte", E, 0.9, ""),
    ("visionkb:ds-rpechoroid-snrna", "Melanocytes"): ("CL:4030000", "choroidal melanocyte", E, 0.9, ""),
}
CURIE_MAP = {
    "CL": "http://purl.obolibrary.org/obo/CL_",
    "skos": "http://www.w3.org/2004/02/skos/core#",
    "semapv": "https://w3id.org/semapv/vocab/",
    "visionkb": "https://purl.visionkb.org/id/",
}
COLS = [
    "subject_id",
    "subject_label",
    "subject_source",
    "annotation_level",
    "predicate_id",
    "object_id",
    "object_label",
    "mapping_justification",
    "confidence",
    "author_id",
    "review_status",
    "comment",
]


def slugify(label: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "-", str(label)).strip("-").lower()


def ols_check(ids: dict[str, str]) -> dict:
    out = {}
    for cid, exp in sorted(ids.items()):
        iri = "http://purl.obolibrary.org/obo/" + cid.replace(":", "_")
        for attempt in range(3):
            try:
                r = requests.get(
                    "https://www.ebi.ac.uk/ols4/api/ontologies/cl/terms", params={"iri": iri}, timeout=30
                ).json()
                break
            except requests.RequestException:
                time.sleep(2 * (attempt + 1))
        terms = r.get("_embedded", {}).get("terms", [])
        got = terms[0]["label"] if terms else None
        obs = terms[0].get("is_obsolete") if terms else None
        out[cid] = {"expected": exp, "ols_label": got, "obsolete": obs, "ok": got == exp and not obs}
    return out


def build(labels: dict) -> list[dict]:
    rows = []
    for ds in sorted(labels):
        slug = ds.removeprefix("visionkb:ds-")
        lv = labels[ds]
        items = [("majorclass", x) for x in lv["majorclass"]] + [
            ("celltype", x) for x in sorted(lv.get("celltype", {}))
        ]
        for level, lab in items:
            m = OVR.get((ds, lab)) or (M if level == "majorclass" else CT).get(lab)
            if m is None:
                continue
            cid, clab, pred, conf, com = m
            rows.append(
                {
                    "subject_id": f"visionkb:cs-{slug}-{slugify(lab)}",
                    "subject_label": lab,
                    "subject_source": ds,
                    "annotation_level": level,
                    "predicate_id": pred,
                    "object_id": cid,
                    "object_label": clab,
                    "mapping_justification": "semapv:ManualMappingCuration",
                    "confidence": conf,
                    "author_id": "visionkb:curators",
                    "review_status": "proposed",
                    "comment": com,
                }
            )
    return rows


def main():
    labels = json.loads(Path(sys.argv[1]).read_text())
    rows = build(labels)
    unmapped_mc = [
        (d, x)
        for d in labels
        for x in labels[d]["majorclass"]
        if not any(
            r["subject_source"] == d and r["annotation_level"] == "majorclass" and r["subject_label"] == x for r in rows
        )
    ]
    if unmapped_mc:
        sys.exit(f"majorclass labels without a mapping: {unmapped_mc}")
    if "--no-ols" not in sys.argv:
        chk = ols_check({r["object_id"]: r["object_label"] for r in rows})
        bad = {k: v for k, v in chk.items() if not v["ok"]}
        if bad:
            sys.exit(f"OLS check failed: {bad}")
        print(f"OLS: {len(chk)} CL ids verified")
    hdr = ["# curie_map:"] + [f"#   {k}: {v}" for k, v in CURIE_MAP.items()]
    hdr += [
        "# mapping_set_id: https://purl.visionkb.org/crosswalk/heca-celltypes/v0.1",
        "# license: https://creativecommons.org/licenses/by/4.0/",
        f"# mapping_date: {datetime.date.today()}",
        "# comment: HECA atlas labels (majorclass and celltype levels) to Cell Ontology; labels without a row are CellStates",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        f.write("\n".join(hdr) + "\n")
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_ct = sum(len(v.get("celltype", {})) for v in labels.values())
    n_ct_mapped = sum(r["annotation_level"] == "celltype" for r in rows)
    print(f"wrote {OUT}: {len(rows)} rows; celltype labels mapped {n_ct_mapped}/{n_ct}")


if __name__ == "__main__":
    main()
