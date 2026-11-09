"""Shared edge construction (contracts §5, §7)."""

from __future__ import annotations

from ..ids import edge_id, qualifiers_of

HUMAN = "NCBITaxon:9606"


def make_edge(
    subject: str, predicate: str, obj: str, *, source: str, knowledge_level: str, agent_type: str, release: str, **props
) -> dict:
    e = {"subject": subject, "predicate": predicate, "object": obj, **props}
    e["id"] = edge_id(subject, predicate, obj, qualifiers_of(e))
    e.update(
        primary_knowledge_source=source,
        knowledge_level=knowledge_level,
        agent_type=agent_type,
        kg_release=release,
    )
    e.setdefault("status", "curated")
    e.setdefault("visibility", "public")
    return e
