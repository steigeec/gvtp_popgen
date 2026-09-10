#!/usr/bin/env python3
# coding: utf-8
# Author : Diana Aguilar
 
import pandas as pd
import argparse
 
#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("fixed_af", type=float, help="frequency of allele to be considered fixed, for low coverage data it is not recommended to use 1")
parser.add_argument("pop1", type=str, help="name of population 1 (prefix of maf file)")
args = parser.parse_args()
 
#Read files
outfile = args.pop1 + ".mafs"
fixed_af = args.fixed_af
pop1=pd.read_csv(args.pop1 + ".mafs", sep="\t", skiprows=[1], dtype={"chromo": str, "position": int, "major": str, "minor": str, "ref": str, "knownEM": float, "nInd": int})
polysites=pd.read_csv("/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog.sites.gz",sep="\t",header=None,names=["chromo","position"], dtype=str)
polysites["position"]=pd.to_numeric(polysites["position"], errors="coerce")
polysites=polysites.dropna().astype({"chromo": str, "position": int})
 
#Merge with confirmed polymorphic sites
pop1poly= pop1.merge(polysites, on=["chromo","position"], how='left',indicator='PolyConfirmed')
pop1poly['PolyBinary'] = pop1poly['PolyConfirmed'].eq('both')
 
#Keep if confirmed polymorphic or fixed
KeepSites=pop1poly[((pop1poly.knownEM<1-fixed_af) | (pop1poly.knownEM>fixed_af)) | (pop1poly.PolyBinary==True)]
KeepSites.to_csv(outfile,columns=["chromo","position","major","minor","ref","knownEM","nInd"],sep="\t",index=False)
