#!/usr/bin/env python3
"""
Build the GO ancestor closure (is_a + part_of only, the true-path relations)
from a current go-basic.obo, used in place of GO.db's ancestor tables, which
lag behind the current ontology.

Outputs (written next to out_path):
  go_ancestor_map.tsv     one row per (term, ancestor) pair
  go_term_names.tsv       id, name, ontology (BP/MF/CC)
  go_slim_plant_terms.txt GO ids tagged `subset: goslim_plant`
  go_obsolete_remap.tsv   old_id, new_id for obsolete terms with a single
                          replaced_by successor (followed transitively);
                          terms with only `consider:` suggestions are not remapped

Usage: python3 build_go_ancestor_map.py [obo_path] [out_path]
  obo_path: path to go-basic.obo (default: downloads fresh copy to /tmp)
  out_path: default BASE_DIR/genome/go_ancestor_map.tsv
"""
import sys
import os
import urllib.request
from collections import defaultdict

BASE_DIR = "/private/groups/meyerlab/esteiger"

OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
DEFAULT_OBO = "/tmp/go-basic.obo"
DEFAULT_OUT = f"{BASE_DIR}/genome/go_ancestor_map.tsv"


NAMESPACE_ABBREV = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}


def parse_obo(path):
    """Returns (id2name, id2ns, parents, replaced_by, slim_plant).
    parents[id] is the set of direct is_a/part_of targets."""
    id2name = {}
    id2ns = {}
    parents = defaultdict(set)
    replaced_by = {}
    slim_plant = set()
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
            elif line.startswith("namespace: ") and cur_id:
                id2ns[cur_id] = NAMESPACE_ABBREV.get(line[len("namespace: "):].strip(), "")
            elif line.startswith("is_a: ") and cur_id:
                target = line[6:].split("!")[0].strip()
                parents[cur_id].add(target)
            elif line.startswith("relationship: part_of ") and cur_id:
                target = line[len("relationship: part_of "):].split("!")[0].strip()
                parents[cur_id].add(target)
            elif line.startswith("replaced_by: ") and cur_id:
                replaced_by[cur_id] = line[len("replaced_by: "):].strip()
            elif line.startswith("subset: goslim_plant") and cur_id:
                slim_plant.add(cur_id)
            elif line.startswith("is_obsolete: true") and cur_id:
                id2name.pop(cur_id, None)
                id2ns.pop(cur_id, None)
                slim_plant.discard(cur_id)
    return id2name, id2ns, parents, replaced_by, slim_plant


def main():
    obo_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OBO
    out_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUT

    if not os.path.exists(obo_path):
        print(f"Downloading {OBO_URL} -> {obo_path}...")
        urllib.request.urlretrieve(OBO_URL, obo_path)

    id2name, id2ns, parents, replaced_by, slim_plant = parse_obo(obo_path)
    print(f"Parsed {len(id2name)} non-obsolete terms, {sum(len(v) for v in parents.values())} direct is_a/part_of edges.")

    # Memoized ancestor closure (GO is a DAG, so plain recursion + memo is safe).
    memo = {}
    import sys as _sys
    _sys.setrecursionlimit(10000)

    def ancestors(term_id, _stack=None):
        if term_id in memo:
            return memo[term_id]
        if _stack is None:
            _stack = set()
        if term_id in _stack:
            return set()  # cycle guard, shouldn't happen in a real DAG
        _stack = _stack | {term_id}
        result = set()
        for p in parents.get(term_id, ()):
            if p not in id2name:
                continue  # obsolete/unknown parent, skip
            result.add(p)
            result |= ancestors(p, _stack)
        memo[term_id] = result
        return result

    n_written = 0
    with open(out_path, "w") as out:
        out.write("term\tancestor\n")
        for term_id in id2name:
            for anc in sorted(ancestors(term_id)):
                out.write(f"{term_id}\t{anc}\n")
                n_written += 1

    print(f"Wrote {n_written} (term, ancestor) pairs for {len(id2name)} terms to {out_path}")

    names_path = os.path.join(os.path.dirname(out_path) or ".", "go_term_names.tsv")
    with open(names_path, "w") as out:
        out.write("id\tname\tontology\n")
        for term_id, name in id2name.items():
            out.write(f"{term_id}\t{name}\t{id2ns.get(term_id, '')}\n")
    print(f"Wrote {len(id2name)} (id, name, ontology) rows to {names_path}")

    slim_path = os.path.join(os.path.dirname(out_path) or ".", "go_slim_plant_terms.txt")
    with open(slim_path, "w") as out:
        for term_id in sorted(slim_plant):
            out.write(f"{term_id}\n")
    print(f"Wrote {len(slim_plant)} goslim_plant term ids to {slim_path}")

    # Follow replaced_by chains until a live term; cycles and dead ends are not remapped.
    def resolve(term_id, _stack=None):
        _stack = _stack or set()
        if term_id in id2name:
            return term_id
        if term_id in _stack or term_id not in replaced_by:
            return None
        return resolve(replaced_by[term_id], _stack | {term_id})

    remap_path = os.path.join(os.path.dirname(out_path) or ".", "go_obsolete_remap.tsv")
    n_remapped = 0
    n_dead_end = 0
    with open(remap_path, "w") as out:
        out.write("old_id\tnew_id\n")
        for old_id in replaced_by:
            new_id = resolve(old_id)
            if new_id is None:
                n_dead_end += 1
                continue
            out.write(f"{old_id}\t{new_id}\n")
            n_remapped += 1
    print(f"Wrote {n_remapped} obsolete->current GO id remaps to {remap_path} "
          f"({n_dead_end} had a replaced_by chain that didn't resolve to a live term)")


if __name__ == "__main__":
    main()
