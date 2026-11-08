#!/usr/bin/env python3
"""
Translate the Deinandra CDS set
(remove_duplicates_rm_internalstops_longest_isoform.cds) to protein for
OrthoFinder. Translation stops at the first in-frame stop codon; records that
still contain internal stops are counted and reported.

Usage: python3 translate_cds_to_pep.py [output.pep]
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq

BASE_DIR = "/private/groups/meyerlab/esteiger"

CDS = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.cds"


def main():
    out_file = sys.argv[1] if len(sys.argv) > 1 else \
        f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.pep"

    n_total = 0
    n_trimmed = 0       # length not a multiple of 3, trailing partial codon dropped
    n_internal_stop = 0  # stop codon found before the last codon
    n_no_stop = 0        # no in-frame stop at all (likely a partial/non-terminated CDS)
    n_empty = 0          # translated to zero residues, excluded from output
    written = []

    for rec in SeqIO.parse(CDS, "fasta"):
        n_total += 1
        seq = str(rec.seq).upper()
        rem = len(seq) % 3
        if rem:
            seq = seq[: len(seq) - rem]
            n_trimmed += 1

        full_aa = str(Seq(seq).translate(table=1, to_stop=False))
        stop_pos = full_aa.find("*")

        if stop_pos == -1:
            protein = full_aa
            n_no_stop += 1
        else:
            protein = full_aa[:stop_pos]
            if stop_pos != len(full_aa) - 1:
                n_internal_stop += 1

        if not protein:
            n_empty += 1
            continue

        written.append((rec.id, protein))

    with open(out_file, "w") as out:
        for gid, protein in written:
            out.write(f">{gid}\n")
            for i in range(0, len(protein), 60):
                out.write(protein[i : i + 60] + "\n")

    print(f"Translated {n_total} CDS records -> {len(written)} protein sequences in {out_file}")
    print(f"  Trailing partial codon trimmed: {n_trimmed}")
    print(f"  Internal (premature) stop codon found: {n_internal_stop}")
    print(f"  No in-frame stop at all (likely non-terminated): {n_no_stop}")
    print(f"  Empty translation (excluded): {n_empty}")


if __name__ == "__main__":
    main()
