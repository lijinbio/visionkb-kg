"""Generate schema/generated/visionkb.schema.json: JSON Schema of the visionkb classes and what they reference.

The full Biolink JSON Schema is ~6 MB; only definitions reachable from visionkb classes are kept.
Usage: python scripts/gen_json_schema.py [--check]
"""

import json
import sys
from pathlib import Path

from linkml.generators.jsonschemagen import JsonSchemaGenerator
from linkml_runtime.utils.schemaview import SchemaView

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "src" / "visionkb_kg" / "schema" / "visionkb.yaml"
OUT = ROOT / "src" / "visionkb_kg" / "schema" / "generated" / "visionkb.schema.json"


def build() -> str:
    sv = SchemaView(str(SCHEMA))
    own = [c.name for c in sv.all_classes(imports=False).values()]
    full = json.loads(JsonSchemaGenerator(str(SCHEMA), mergeimports=True).serialize())
    defs = full.get("$defs", {})
    keep, todo = set(), [n.replace(" ", "") for n in own]
    names = {k.lower(): k for k in defs}
    todo = [names.get(n.lower(), n) for n in todo]
    while todo:
        n = todo.pop()
        if n in keep or n not in defs:
            continue
        keep.add(n)
        stack = [defs[n]]
        while stack:
            x = stack.pop()
            if isinstance(x, dict):
                ref = x.get("$ref")
                if isinstance(ref, str) and ref.startswith("#/$defs/"):
                    todo.append(ref.split("/")[-1])
                stack.extend(x.values())
            elif isinstance(x, list):
                stack.extend(x)
    out = {k: v for k, v in full.items() if k not in ("$defs", "properties", "additionalProperties", "required")}
    out["title"] = "visionkb"
    out["description"] = (
        f"VisionKB classes ({', '.join(sorted(own))}) on Biolink 4.4.6; definitions pruned to those reachable."
    )
    out["$defs"] = {k: defs[k] for k in sorted(keep)}
    # linkml emits the multivalued type designator `category` as {type: array, enum: [...]}, which no array can satisfy;
    # require instead that the list contains the class CURIE (KGX lists ancestors in category too).
    for d in out["$defs"].values():
        cat = d.get("properties", {}).get("category")
        if isinstance(cat, dict) and cat.get("type") == "array" and "enum" in cat:
            cat["contains"] = {"enum": cat.pop("enum")}
    return json.dumps(out, indent=2, sort_keys=True) + "\n"


if __name__ == "__main__":
    text = build()
    if "--check" in sys.argv:
        if OUT.read_text() != text:
            sys.exit(f"{OUT} is out of date: run python scripts/gen_json_schema.py")
    else:
        OUT.write_text(text)
        print(f"wrote {OUT} ({len(text):,} bytes, {len(json.loads(text)['$defs'])} definitions)")
