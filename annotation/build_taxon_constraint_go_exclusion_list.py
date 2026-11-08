#!/usr/bin/env python3
"""
Build a GO term exclusion list from the GO Consortium's never_in_taxon
Viridiplantae (NCBITaxon:33090) assertions, which are present only in
go-plus.owl (RO_0002161), not in go.obo.
  1. Stream go-plus.owl and collect GO terms asserted never_in_taxon Viridiplantae.
  2. Add every is_a/part_of descendant of those terms, using go.obo.
only_in_taxon constraints are not handled.

Usage: python3 build_taxon_constraint_go_exclusion_list.py [owl_path] [obo_path] [out_path]
  owl_path: path to go-plus.owl (default: downloads fresh copy to /tmp)
  obo_path: path to go.obo, is_a/part_of only (default: downloads fresh copy to /tmp)
  out_path: default BASE_DIR/genome/go_terms_excluded_taxon_constraint.txt
"""
import sys
import os
import urllib.request
from collections import defaultdict
import xml.etree.ElementTree as ET

BASE_DIR = "/private/groups/meyerlab/esteiger"

OWL_URL = "https://current.geneontology.org/ontology/extensions/go-plus.owl"
OBO_URL = "https://current.geneontology.org/ontology/go.obo"
DEFAULT_OWL = "/tmp/go-plus.owl"
DEFAULT_OBO = "/tmp/go_full.obo"
DEFAULT_OUT = f"{BASE_DIR}/genome/go_terms_excluded_taxon_constraint.txt"

RDF_NS = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
OWL_NS = "{http://www.w3.org/2002/07/owl#}"
OBO_NS = "{http://purl.obolibrary.org/obo/}"

NEVER_IN_TAXON_TAG = OBO_NS + "RO_0002161"
VIRIDIPLANTAE_URI = "http://purl.obolibrary.org/obo/NCBITaxon_33090"


def uri_to_curie(uri):
    """http://purl.obolibrary.org/obo/GO_0007626 -> GO:0007626"""
    frag = uri.rsplit("/", 1)[-1]
    return frag.replace("_", ":", 1)


def extract_never_in_plants(owl_path):
    """Return GO curies directly asserted never_in_taxon Viridiplantae."""
    seeds = set()
    n_classes = 0
    context = ET.iterparse(owl_path, events=("end",))
    for _, elem in context:
        if elem.tag != OWL_NS + "Class":
            continue
        n_classes += 1
        about = elem.get(RDF_NS + "about", "")
        if "/GO_" in about:
            for child in elem.findall(NEVER_IN_TAXON_TAG):
                if child.get(RDF_NS + "resource") == VIRIDIPLANTAE_URI:
                    seeds.add(uri_to_curie(about))
                    break
        elem.clear()
    print(f"Scanned {n_classes} owl:Class elements; "
          f"found {len(seeds)} GO terms directly asserted never_in_taxon Viridiplantae")
    return seeds


def parse_obo_is_a_part_of(path):
    """Returns (id2name, parents) using is_a/part_of edges only; a taxon
    constraint on X does not carry over to "regulation of X"."""
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
            elif line.startswith("is_obsolete: true") and cur_id:
                id2name.pop(cur_id, None)
    return id2name, parents


def main():
    owl_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OWL
    obo_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OBO
    out_path = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_OUT

    if not os.path.exists(owl_path):
        print(f"Downloading {OWL_URL} -> {owl_path}...")
        urllib.request.urlretrieve(OWL_URL, owl_path)
    if not os.path.exists(obo_path):
        print(f"Downloading {OBO_URL} -> {obo_path}...")
        urllib.request.urlretrieve(OBO_URL, obo_path)

    seeds = extract_never_in_plants(owl_path)
    id2name, parents = parse_obo_is_a_part_of(obo_path)

    children = defaultdict(set)
    for child, pars in parents.items():
        for p in pars:
            children[p].add(child)

    # seeds absent from the obo graph are still excluded themselves
    excluded = set()
    frontier = set(seeds)
    while frontier:
        excluded |= frontier
        next_frontier = set()
        for t in frontier:
            next_frontier |= children.get(t, set())
        frontier = next_frontier - excluded

    with open(out_path, "w") as f:
        for gid in sorted(excluded):
            f.write(f"{gid}\t{id2name.get(gid, '')}\n")

    print(f"{len(seeds)} directly-asserted seed terms -> {len(excluded)} total "
          f"(including is_a/part_of descendants)")
    print(f"Wrote {len(excluded)} taxon-constraint-excluded (never-in-Viridiplantae) "
          f"GO terms to {out_path}")


if __name__ == "__main__":
    main()
