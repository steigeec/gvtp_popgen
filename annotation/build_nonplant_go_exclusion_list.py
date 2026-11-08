#!/usr/bin/env python3
"""
Build a list of GO term IDs implausible for a plant (mostly animal-specific
terms carried over by BLAST-based GO transfer), to filter from GO enrichment
results. Excludes every descendant of a few high-level animal-only anchor
terms, walking is_a/part_of/regulates edges parsed from a current
go-basic.obo. A heuristic filter, not formal taxon-constraint curation.

Usage: python3 build_nonplant_go_exclusion_list.py [obo_path] [out_path]
  obo_path: path to go-basic.obo (default: downloads fresh copy to /tmp)
  out_path: default BASE_DIR/genome/go_terms_excluded_nonplant.txt
"""
import sys
import re
import urllib.request
from collections import defaultdict

BASE_DIR = "/private/groups/meyerlab/esteiger"

OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
DEFAULT_OBO = "/tmp/go-basic.obo"
DEFAULT_OUT = f"{BASE_DIR}/genome/go_terms_excluded_nonplant.txt"

# High-level terms with no plausible plant equivalent
ANCHOR_NAMES = {
    "animal organ morphogenesis",
    "animal organ regeneration",
    "animal organ development",
    "animal organ formation",
    "animal organ maturation",
    "animal gross anatomical part developmental process",
    "nervous system development",
    "nervous system process",
    "synapse",
    "synaptic signaling",       # not a descendant of "synapse", needs its own anchor
    "immune system process",
    "circulatory system process",
    "muscle structure development",
    "muscle system process",
    "contractile muscle fiber",
    "sensory perception",
    "digestion",
    "phagocytosis",   # also carried by Viridiplantae-scoped eggNOG annotations
    "transcytosis",
}

# Traversal stops at these nodes (neither they nor their descendants are
# excluded): they sit under an animal anchor but their descendants include
# plant terms (flower/sepal development, phloem transport, plant embryogenesis).
BLOCK_NAMES = {
    "reproductive system development",
    "reproductive structure development",
    "vascular transport",
    "embryo development",
}


REGULATES_PREFIXES = (
    "relationship: regulates ",
    "relationship: positively_regulates ",
    "relationship: negatively_regulates ",
)


def parse_obo(path):
    """Returns (id2name, parents) where parents[id] = set of direct
    is_a/part_of/regulates targets. Unlike build_go_ancestor_map.py, regulates
    edges are followed so that "regulation of X" terms are excluded with X."""
    id2name = {}
    parents = defaultdict(set)
    cur_id = None
    in_term = False
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                in_term = True
                cur_id = None
                continue
            if line.startswith("[") and line != "[Term]":
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                cur_id = line[4:].strip()
            elif line.startswith("name: ") and cur_id:
                id2name[cur_id] = line[6:].strip()
            elif line.startswith("is_a: ") and cur_id:
                target = line[6:].split("!")[0].strip()
                parents[cur_id].add(target)
            elif line.startswith("relationship: part_of ") and cur_id:
                target = line[len("relationship: part_of "):].split("!")[0].strip()
                parents[cur_id].add(target)
            elif cur_id and line.startswith(REGULATES_PREFIXES):
                prefix = next(p for p in REGULATES_PREFIXES if line.startswith(p))
                target = line[len(prefix):].split("!")[0].strip()
                parents[cur_id].add(target)
            elif line.startswith("is_obsolete: true") and cur_id:
                id2name.pop(cur_id, None)
    return id2name, parents


def main():
    obo_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OBO
    out_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUT

    import os
    if not os.path.exists(obo_path):
        print(f"Downloading {OBO_URL} -> {obo_path}...")
        urllib.request.urlretrieve(OBO_URL, obo_path)

    id2name, parents = parse_obo(obo_path)
    name2id = {v: k for k, v in id2name.items()}

    anchors = set()
    for name in ANCHOR_NAMES:
        gid = name2id.get(name)
        if gid is None:
            print(f"WARNING: anchor term not found in ontology: {name!r}", file=sys.stderr)
            continue
        anchors.add(gid)
    print(f"Resolved {len(anchors)}/{len(ANCHOR_NAMES)} anchor terms.")

    blocked = set()
    for name in BLOCK_NAMES:
        gid = name2id.get(name)
        if gid is None:
            print(f"WARNING: block term not found in ontology: {name!r}", file=sys.stderr)
            continue
        blocked.add(gid)
    print(f"Resolved {len(blocked)}/{len(BLOCK_NAMES)} block terms.")

    # BFS from anchors over the reversed graph
    children = defaultdict(set)
    for child, pars in parents.items():
        for p in pars:
            children[p].add(child)

    excluded = set()
    frontier = set(anchors) - blocked
    while frontier:
        excluded |= frontier
        next_frontier = set()
        for t in frontier:
            if t in blocked:
                continue  # don't expand past a blocked node
            next_frontier |= children.get(t, set())
        frontier = next_frontier - excluded - blocked

    with open(out_path, "w") as f:
        for gid in sorted(excluded):
            f.write(f"{gid}\t{id2name.get(gid, '')}\n")

    print(f"Wrote {len(excluded)} excluded (non-plant-plausible) GO terms to {out_path}")


if __name__ == "__main__":
    main()
