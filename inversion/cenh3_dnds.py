#!/usr/bin/env python3
"""
dN/dS of CENH3 in D. increscens subsp. villosa versus relatives, and against a
genome-wide background of Deinandra-sunflower 1:1 orthologs.

CENH3 sequences (coding sequences; cds/ and OUT_DIR/cenh3_cds_raw.fa):
  Deinandra       IDmodified-mRNA-22870, a gene model that merges an Argonaute
                  gene with CENH3; the CENH3 part (DEINVI_00019:23,147,975-
                  23,151,333, + strand) is the model's last 175 codons, trimmed
                  here to start at its Met (MARTKHP...). The only CENH3 in the
                  assembly (tblastn of sunflower CENH3 tails finds no other).
  H. annuus       HanXRQr2_Chr14g0649211, HanXRQr2_Chr15g0703621 (Ensembl Plants)
  H. anomalus     MFS8009295.1 (Chr14), MFS8018859.1 (Chr15) (NCBI)
  H. debilis      KAM0005100.1 (Chr14-type; NCBI, GCA_052426805.1); no Chr15-type copy found
  Layia glandulosa  miniprot gene model on the alternate-haplotype assembly
                  GCA_046128875.1 (JBJXVY010000114.1:879,022-882,030; the primary
                  assembly lacks CENH3); a near-identical second copy is not used
  (sequence set assembled by cenh3_collect_cds.py)
  L. sativa       KAJ0221726.1 (chr2): histone-fold domain only -- its annotated
                  N-terminus is not a CENH3 tail (mispredicted start), so
                  codons before the histone fold are masked. KAJ0205668.1 (chr5)
                  is fragmentary (probable pseudogene) and excluded.
  A. thaliana     AT1G01370.1 (outgroup)
  Tomato was not retrievable (Ensembl Plants CDS unavailable at run time).

Steps
  1. Protein alignment (MAFFT L-INS-i), codon alignment (PAL2NAL).
  2. Domains: N-terminal tail = alignment columns before the conserved
     histone-fold start (the R of the R[FY][KR]PGT motif, located in H. annuus Chr14
     because Deinandra carries RFKAGT);
     histone-fold domain (HFD) = from there to the end.
  3. Pairwise dN, dS, dN/dS (yn00, Yang & Nielsen 2000) for every pair, for the
     full protein, tail and HFD.
  4. codeml site models on the full alignment: M1a vs M2a and M7 vs M8
     (likelihood-ratio tests, chi-square 2 df), with Bayes Empirical Bayes
     sites under M8; branch-site model A on the Deinandra branch vs its null
     (chi-square 1 df; conservative 50:50 mixture p also reported).
     Tree: (Athaliana, Lsativa_chr2, ((Deinandra, Layia), ((Hannuus_Chr14,
     Hanomalus_Chr14, Hdebilis_Chr14), (Hannuus_Chr15, Hanomalus_Chr15)))) --
     Deinandra and Layia are both Madiinae; the two
     Helianthus copies pair by chromosome across species; placing the
     Helianthus duplication after the split from Deinandra is checked in the
     yn00 output (paralog dS vs Deinandra-Helianthus dS).
  5. Background: N_BACKGROUND random Deinandra-H. annuus 1:1 orthologs
     (GENESPACE syntenic anchors in the same orthogroup), each aligned and run
     through yn00 the same way; CENH3's Deinandra-H. annuus dN/dS is placed in
     that distribution (percentile).

Usage: cenh3_dnds.py [n_procs]
Output: OUT_DIR/ (alignments, yn00_pairs_cenh3.tsv, codeml_tests.tsv,
        codeml_M8_BEB_sites.tsv, background_yn00.tsv, summary.txt)
Conda env: dnds (BASE_DIR/software/miniconda3/envs/dnds)
  conda create -n dnds -c conda-forge -c bioconda paml=4.10.10 mafft=7.526 \
    pal2nal=14.1 miniprot=0.18
(the tools are called by path through DNDS_ENV)
"""

import gzip
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from scipy.stats import chi2

BASE_DIR = "/private/groups/meyerlab/esteiger"

OUT_DIR = f"{BASE_DIR}/inversions/centromere/cenh3_dnds"
CDS_DIR = f"{OUT_DIR}/cds"
DNDS_ENV = f"{BASE_DIR}/software/miniconda3/envs/dnds/bin"
MAFFT, PAL2NAL, YN00, CODEML = (f"{DNDS_ENV}/{t}" for t in ("mafft", "pal2nal.pl", "yn00", "codeml"))
DEIN_CDS = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.cds"
HA_CDS = f"{CDS_DIR}/Helianthus_annuus.HanXRQr2.0-SUNRISE.cds.all.fa.gz"
SYNHITS = f"{BASE_DIR}/inversions/synteny/genespace/syntenicHits/Helianthus_annuus_vs_Deinandra.synHits.txt.gz"
N_BACKGROUND, SEED = 2000, 1
TREE = ("(Athaliana,Lsativa_chr2,((Deinandra,Layia),"
        "((Hannuus_Chr14,Hanomalus_Chr14,Hdebilis_Chr14),(Hannuus_Chr15,Hanomalus_Chr15))));")
TREE_FG = TREE.replace("Deinandra,", "Deinandra #1,")


def read_fasta(path, keep=lambda h: True):
    op = gzip.open if path.endswith(".gz") else open
    out, cur = {}, None
    with op(path, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                cur = l[1:].strip() if keep(l[1:].strip()) else None
                if cur:
                    out[cur] = []
            elif cur:
                out[cur].append(l.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def clean_cds(s):
    s = s[:len(s) - len(s) % 3]
    if str(Seq(s[-3:]).translate()) == "*":
        s = s[:-3]
    return s


def run(cmd, cwd=None, stdout=None, check=True):
    subprocess.run(cmd, cwd=cwd, stdout=stdout or subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=check)


def codon_align(cds, workdir, linsi=True):
    """cds: {name: nt}; returns codon alignment {name: aligned nt} (PAL2NAL, gaps kept)"""
    with open(f"{workdir}/p.fa", "w") as f:
        f.writelines(f">{k}\n{Seq(v).translate()}\n" for k, v in cds.items())
    with open(f"{workdir}/n.fa", "w") as f:
        f.writelines(f">{k}\n{v}\n" for k, v in cds.items())
    with open(f"{workdir}/p.aln", "w") as fh:
        run([MAFFT] + (["--localpair", "--maxiterate", "1000"] if linsi else ["--auto"]) + ["--quiet", f"{workdir}/p.fa"], stdout=fh)
    with open(f"{workdir}/c.aln", "w") as fh:
        run([PAL2NAL, f"{workdir}/p.aln", f"{workdir}/n.fa", "-output", "fasta"], stdout=fh)
    return read_fasta(f"{workdir}/c.aln")


def write_paml(aln, path):
    L = len(next(iter(aln.values())))
    with open(path, "w") as f:
        f.write(f" {len(aln)} {L}\n")
        for k, v in aln.items():
            f.write(f"{k}\n{v}\n")


def yn00(aln, workdir):
    """Pairwise dN, dS, omega (Yang & Nielsen 2000) as a DataFrame."""
    write_paml(aln, f"{workdir}/seq.phy")
    with open(f"{workdir}/yn00.ctl", "w") as f:
        f.write("seqfile = seq.phy\noutfile = yn.out\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n")
    run([YN00, "yn00.ctl"], cwd=workdir)
    names = list(aln)
    txt = open(f"{workdir}/yn.out").read()
    block = txt.split("Yang & Nielsen (2000) method")[1].split("LWL85")[0] if "LWL85" in txt else \
        txt.split("Yang & Nielsen (2000) method")[1]
    rows = []
    for l in block.splitlines():
        p = l.split()
        if len(p) >= 11 and p[0].isdigit() and p[1].isdigit():
            i, j = int(p[0]) - 1, int(p[1]) - 1
            rows.append(dict(seq1=names[i], seq2=names[j], S=float(p[2]), N=float(p[3]), omega=float(p[6]),
                             dN=float(p[7]), dN_SE=float(p[9]), dS=float(p[10]), dS_SE=float(p[12])))
    return pd.DataFrame(rows)


def codeml(aln, workdir, model, tree):
    """model: M1a, M2a, M7, M8, A (branch-site), Anull; returns (lnL, np, mlc text)"""
    settings = {"M1a": "model = 0\nNSsites = 1\nfix_omega = 0\nomega = 0.5",
                "M2a": "model = 0\nNSsites = 2\nfix_omega = 0\nomega = 0.5",
                "M7": "model = 0\nNSsites = 7\nfix_omega = 0\nomega = 0.5",
                "M8": "model = 0\nNSsites = 8\nfix_omega = 0\nomega = 0.5",
                "A": "model = 2\nNSsites = 2\nfix_omega = 0\nomega = 1.5",
                "Anull": "model = 2\nNSsites = 2\nfix_omega = 1\nomega = 1"}[model]
    d = f"{workdir}/{model}"
    os.makedirs(d, exist_ok=True)
    write_paml(aln, f"{d}/seq.phy")
    open(f"{d}/tree.nwk", "w").write(tree + "\n")
    open(f"{d}/codeml.ctl", "w").write(
        f"seqfile = seq.phy\ntreefile = tree.nwk\noutfile = mlc\nnoisy = 0\nverbose = 0\nrunmode = 0\n"
        f"seqtype = 1\nCodonFreq = 2\nclock = 0\nfix_kappa = 0\nkappa = 2\n{settings}\n"
        f"cleandata = 0\nfix_blength = 0\nmethod = 0\n")
    # this codeml build exits 1 even on success; check for the lnL line instead
    run([CODEML, "codeml.ctl"], cwd=d, check=False)
    mlc = open(f"{d}/mlc").read()
    if "lnL(ntime" not in mlc:
        raise RuntimeError(f"codeml {model} produced no likelihood (see {d}/mlc)")
    m = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?[\d.]+)", mlc)
    return float(m.group(2)), int(m.group(1)), mlc


def beb_sites(mlc, aln_ref):
    """Bayes Empirical Bayes sites (M8) with P(omega > 1) > 0.5, with residue in the reference sequence"""
    if "Bayes Empirical Bayes" not in mlc:
        return pd.DataFrame()
    sec = mlc.split("Bayes Empirical Bayes")[1].split("Positively selected sites")[1].split("The grid")[0]
    rows = []
    for l in sec.splitlines():
        m = re.match(r"\s*(\d+)\s+(\S)\s+([\d.]+)(\*{0,2})\s+([\d.]+)\s+\+-\s+([\d.]+)", l)
        if m:
            rows.append(dict(codon_column=int(m.group(1)), residue_first_seq=m.group(2),
                             P_omega_gt_1=float(m.group(3)), signif=m.group(4), post_mean_omega=float(m.group(5))))
    return pd.DataFrame(rows)


def background_pair(args):
    """yn00 for one Deinandra - H. annuus ortholog pair"""
    gid, d_nt, h_nt = args
    with tempfile.TemporaryDirectory(dir="/data/tmp") as w:
        try:
            aln = codon_align({"D": d_nt, "H": h_nt}, w, linsi=False)
            r = yn00(aln, w)
            if len(r):
                return dict(pair=gid, length_codons=len(d_nt) // 3, **r.iloc[0][["S", "N", "omega", "dN", "dS"]].to_dict())
        except Exception:
            return None


def main():
    n_procs = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    os.makedirs(OUT_DIR, exist_ok=True)
    work = tempfile.mkdtemp(dir="/data/tmp")

    # ---- 1. CENH3 set ----
    raw = read_fasta(f"{OUT_DIR}/cenh3_cds_raw.fa")
    cds = {}
    for k, v in raw.items():
        if k == "Lsativa_chr5":
            continue
        v = clean_cds(v)
        p = str(Seq(v).translate())
        if k == "Deinandra":
            v = v[3 * p.rindex("MARTK"):]           # last Met-start: the CENH3 part of the merged model
        if k == "Lsativa_chr2":                      # keep histone-fold part only
            v = v[3 * p.index("KPHRF"):]
        cds[k] = v
    aln = codon_align(cds, work, linsi=True)
    with open(f"{OUT_DIR}/cenh3_codon_aln.fa", "w") as f:
        f.writelines(f">{k}\n{v}\n" for k, v in aln.items())
    prot = {k: str(Seq(v.replace("-", "N")).translate()).replace("X", "-") for k, v in aln.items()}
    with open(f"{OUT_DIR}/cenh3_protein_aln.fa", "w") as f:
        f.writelines(f">{k}\n{v}\n" for k, v in prot.items())

    # ---- 2. domains ----
    dp = prot["Deinandra"]
    # boundary located in H. annuus Chr14 (canonical RFKPGT; Deinandra has RFKAGT)
    m = re.search(r"R-*[FY]-*[KR]-*[PA]-*G-*T", prot["Hannuus_Chr14"])
    hfd_col = m.start()
    tail = {k: v[:3 * hfd_col] for k, v in aln.items() if k != "Lsativa_chr2"}
    hfd = {k: v[3 * hfd_col:] for k, v in aln.items()}
    print(f"alignment: {len(dp)} codon columns; N-terminal tail = columns 1-{hfd_col}, histone fold = {hfd_col + 1}-{len(dp)}")
    for k, v in prot.items():
        print(f"  {k:16s} {v[:hfd_col]} | {v[hfd_col:]}")

    # ---- 3. pairwise yn00 ----
    res = []
    for region, a in [("full", aln), ("N-terminal tail", tail), ("histone fold", hfd)]:
        w = tempfile.mkdtemp(dir=work)
        r = yn00(a, w); r.insert(0, "region", region); res.append(r)
    yn = pd.concat(res)
    yn.to_csv(f"{OUT_DIR}/yn00_pairs_cenh3.tsv", sep="\t", index=False)

    # ---- 4. codeml ----
    tests = []
    fits = {mdl: codeml(aln, work, mdl, TREE_FG if mdl.startswith("A") else TREE) for mdl in ("M1a", "M2a", "M7", "M8", "Anull", "A")}
    for h0, h1, df in [("M1a", "M2a", 2), ("M7", "M8", 2), ("Anull", "A", 1)]:
        lr = max(0.0, 2 * (fits[h1][0] - fits[h0][0]))
        p = chi2.sf(lr, df)
        tests.append(dict(null=h0, alternative=h1, lnL_null=fits[h0][0], lnL_alt=fits[h1][0], LRT=lr, df=df,
                          p=p, p_mixture=(p / 2 if h0 == "Anull" else np.nan)))
    tests = pd.DataFrame(tests)
    tests.to_csv(f"{OUT_DIR}/codeml_tests.tsv", sep="\t", index=False)
    beb = beb_sites(fits["M8"][2], aln)
    if len(beb):
        beb["Deinandra_residue"] = [dp[c - 1] for c in beb.codon_column]
        beb["domain"] = np.where(beb.codon_column <= hfd_col, "N-terminal tail", "histone fold")
    beb.to_csv(f"{OUT_DIR}/codeml_M8_BEB_sites.tsv", sep="\t", index=False)
    om8 = re.search(r"p=\s*([\d.]+)\s+q=\s*([\d.]+)\s+\(p1=\s*([\d.]+)\)\s+w=\s*([\d.]+)", fits["M8"][2])

    # ---- 5. background ----
    syn = pd.read_csv(SYNHITS, sep="\t", low_memory=False, usecols=["id1", "id2", "isAnchor", "sameOG", "chr2"])
    syn = syn[syn.isAnchor & syn.sameOG].drop_duplicates("id1").drop_duplicates("id2")
    random.seed(SEED)
    pick = syn.sample(min(N_BACKGROUND * 2, len(syn)), random_state=SEED)
    dcds = read_fasta(DEIN_CDS, lambda h: h.split()[0] in set(pick.id2))
    dcds = {k.split()[0]: v for k, v in dcds.items()}
    hcds = read_fasta(HA_CDS, lambda h: re.search(r"gene:(\S+)", h).group(1) in set(pick.id1))
    hcds = {re.search(r"gene:(\S+)", k).group(1): v for k, v in hcds.items()}
    jobs = []
    for _, r in pick.iterrows():
        if r.id2 in dcds and r.id1 in hcds:
            d, h = clean_cds(dcds[r.id2]), clean_cds(hcds[r.id1])
            if len(d) >= 150 and len(h) >= 150 and "*" not in str(Seq(d).translate()) and "*" not in str(Seq(h).translate()):
                jobs.append((f"{r.id2}|{r.id1}", d, h))
        if len(jobs) >= N_BACKGROUND:
            break
    with Pool(n_procs) as pool:
        bg = pd.DataFrame([x for x in pool.map(background_pair, jobs) if x])
    bg.to_csv(f"{OUT_DIR}/background_yn00.tsv", sep="\t", index=False)
    ok = bg[(bg.dS > 0.01) & (bg.dS < 2) & (bg.omega < 10)]

    # ---- summary ----
    lines = [f"CENH3 dN/dS summary ({len(aln)} sequences; tail columns 1-{hfd_col}, histone fold {hfd_col + 1}-{len(dp)})", ""]
    lines.append("Pairwise yn00, Deinandra vs each relative:")
    sel = yn[(yn.seq1 == "Deinandra") | (yn.seq2 == "Deinandra")]
    for _, r in sel.iterrows():
        other = r.seq2 if r.seq1 == "Deinandra" else r.seq1
        lines.append(f"  {r.region:16s} vs {other:16s} dN {r.dN:.4f}  dS {r.dS:.4f}  dN/dS {r.omega:.3f}")
    par = yn[(yn.region == "full") & (((yn.seq1 == "Hannuus_Chr14") & (yn.seq2 == "Hannuus_Chr15")) | ((yn.seq1 == "Hannuus_Chr15") & (yn.seq2 == "Hannuus_Chr14")))]
    if len(par):
        lines.append(f"  H. annuus paralogs Chr14 vs Chr15 (full): dS {par.dS.iloc[0]:.4f} (tree check: should be < Deinandra-Helianthus dS)")
    lines += ["", "codeml likelihood-ratio tests:"]
    for _, r in tests.iterrows():
        extra = f" (50:50 mixture p = {r.p_mixture:.3g})" if r.null == "Anull" else ""
        lines.append(f"  {r.null} vs {r.alternative}: 2dlnL = {r.LRT:.2f}, df {r.df}, p = {r.p:.3g}{extra}")
    if om8:
        lines.append(f"  M8: proportion of sites in omega>1 class p1 = {om8.group(3)}, omega = {om8.group(4)}")
    if len(beb):
        lines.append(f"  M8 BEB sites with P(omega>1) > 0.5: {len(beb)} ({(beb.P_omega_gt_1 > 0.95).sum()} with P > 0.95)")
    ce = sel[(sel.region == "full") & ((sel.seq1.str.startswith("Hannuus")) | (sel.seq2.str.startswith("Hannuus")))]
    lines += ["", f"Background: {len(ok)} Deinandra-H. annuus 1:1 orthologs (0.01 < dS < 2): "
                  f"median dN/dS {ok.omega.median():.3f} (IQR {ok.omega.quantile(.25):.3f}-{ok.omega.quantile(.75):.3f}), median dS {ok.dS.median():.3f}"]
    for _, r in ce.iterrows():
        other = r.seq2 if r.seq1 == "Deinandra" else r.seq1
        lines.append(f"  CENH3 Deinandra vs {other}: dN/dS {r.omega:.3f} -> percentile {(ok.omega < r.omega).mean() * 100:.1f}; "
                     f"dS {r.dS:.3f} (background percentile {(ok.dS < r.dS).mean() * 100:.1f})")
    open(f"{OUT_DIR}/summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
