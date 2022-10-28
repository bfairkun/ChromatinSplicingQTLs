#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Add_UnstandardizedBetaCol
# @created     : Thursday Oct 27, 2022 11:44:02 CDT
#
# @description : 
######################################################################

import sys
import pysam
import pandas as pd
from Bio import motifs
from Bio.Seq import Seq

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below

def GetRowFromNominalPassTabix(TabixFileObj, var_chrom, var_start, var_id, phe_id):
    fetch = TabixFileObj.fetch(var_chrom, var_start-1, var_start+1)
    for row in fetch:
        if (row[0] == phe_id) & (row[7]==var_id):
            return row

def GetBetaAndP(*args):
    try:
        row = GetRowFromNominalPassTabix(*args)
        return row[11], row[13]
    except:
        return None, None


def GetSeqFromRow(row, RefGenome=None, var_id=None):
    if row['phe_strd'] == '+':
        return str(Seq(RefGenome.fetch(row['phe_chr'], row['phe_from'], row['phe_to']+1)))
    elif row['phe_strd'] == '-':
        return str(Seq(RefGenome.fetch(row['phe_chr'], row['phe_from']-2, row['phe_to']-1)).reverse_complement())
NominalPass_df.apply(GetSeqFromRow, axis=1, RefGenome=RefGenome)

def GetIntersectingGenes(row, TabixFile):
    return [ entry[3] for entry in TabixFile.fetch(row["phe_chr"], row["phe_from"], row["phe_to"], parser=pysam.asTuple())]

def GetVarFromRow(row, RefGenome=None, SwapFromVarID=False):
    if row['phe_strd'] == '+':
        Coords = row['phe_from'], row['phe_to']+1
    elif row['phe_strd'] == '-':
        Coords = row['phe_from']-2, row['phe_to']-1
    if SwapFromVarID:
        try:
            chrom,pos,ref,alt= row['var_id'].split(':')
            # print(Coords, ref, RefGenome.fetch(row['var_chr'], row['var_from'], row['var_to']))
            IsSNP = len(ref) == 1 and len(alt) == 1
            IsInBounds = (row['var_from'] >= Coords[0]) and (row['var_to'] <= Coords[1]) and (row['phe_chr'] == row['var_chr'])
            if not (IsSNP and IsInBounds):
                return None
            LeftSeq = RefGenome.fetch(row['var_chr'], Coords[0], row['var_from']-1)
            VarSeq = RefGenome.fetch(row['var_chr'], row['var_from']-1, row['var_to'])
            # print(VarSeq)
            VarSeq = alt
            RightSeq = RefGenome.fetch(row['var_chr'], row['var_to'], Coords[1])
        except:
            return None
        if row['phe_strd'] == '+':
            RefSeq = LeftSeq + VarSeq + RightSeq
        elif row['phe_strd'] == '-':
            RefSeq = Seq(LeftSeq + VarSeq + RightSeq).reverse_complement()
    else:
        if row['phe_strd'] == '+':
            RefSeq = Seq(RefGenome.fetch(row['phe_chr'], *Coords ))
        elif row['phe_strd'] == '-':
            RefSeq = Seq(RefGenome.fetch(row['phe_chr'], *Coords)).reverse_complement()
    return str(RefSeq)
# NominalPass_df.head()[['var_id', 'phe_strd','phe_from','phe_to', 'var_from', 'var_to']]
# NominalPass_df.head().apply(GetVarFromRow, axis=1, RefGenome=RefGenome, SwapFromVarID=True)
# NominalPass_df.head().apply(GetSeqFromRow, axis=1, RefGenome=RefGenome)





# Create 5'ss motif PWM
Donors = pd.read_csv("SplicingAnalysis/regtools_annotate_combined/basic.5ss.bed.gz", sep='\t', names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
Donors = Donors.loc[Donors['name']=='SpliceDonor_1']
Donors['seq'] = Donors['score'].str[-9:]
m  = motifs.create(Donors['seq'])
pwm = m.counts.normalize(pseudocounts=1)
pssm = pwm.log_odds()

# Get ref genome
RefGenome = pysam.FastaFile("ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa")

files = {
        "polyA_5ss":"QTLs/QTLTools/polyA.Splicing.5PrimeSS.Subset_YRI/NominalPass.txt.gz",
        "chRNA_5ss":"QTLs/QTLTools/chRNA.Splicing.5PrimeSS/NominalPass.txt.gz"
        }

NominalPass_df = pd.concat([pd.read_csv(v, sep=' ').assign(New=k) for k,v in files.items()])

# Add 5'ss seq columns
NominalPass_df['RefSeq'] = NominalPass_df.apply(GetVarFromRow, axis=1, RefGenome=RefGenome, SwapFromVarID=False)
NominalPass_df['AltSeq'] = NominalPass_df.apply(GetVarFromRow, axis=1, RefGenome=RefGenome, SwapFromVarID=True)

# Add 5'ss score columns
NominalPass_df['RefScore'] = NominalPass_df['RefSeq'].apply(pssm.calculate)
NominalPass_df['AltScore'] = NominalPass_df['AltSeq'].apply(lambda x: pssm.calculate(x) if x else None)
NominalPass_df['DeltaPWM'] = NominalPass_df['AltScore'] - NominalPass_df['RefScore']

# Add columns for host-gene
genes_tbx = pysam.TabixFile("ExpressionAnalysis/polyA/ExpressedGeneList.bed.gz")
NominalPass_df["HostGenes"] = NominalPass_df.apply(GetIntersectingGenes, axis=1, TabixFile=genes_tbx)
NominalPass_df = NominalPass_df.explode("HostGenes")

# Add columns for eQTL effect
polyA_eQTL_tabix = pysam.Tabixfile("QTLs/QTLTools/Expression.Splicing.Subset_YRI/NominalPassForColoc.txt.tabix.gz", parser=pysam.asTuple())
NominalPass_df[["polyA_eQTL_P","polyA_eQTL_beta"]] = NominalPass_df.apply(lambda x: GetBetaAndP(polyA_eQTL_tabix, x['var_chr'], x['var_from'], x['var_id'], f"{x['HostGenes']}:{x['HostGenes']}"), axis=1, result_type="expand")

chRNA_eQTL_tabix = pysam.Tabixfile("QTLs/QTLTools/chRNA.Expression.Splicing/NominalPassForColoc.txt.tabix.gz", parser=pysam.asTuple())
NominalPass_df[["chRNA_eQTL_P","chRNA_eQTL_beta"]] = NominalPass_df.apply(lambda x: GetBetaAndP(chRNA_eQTL_tabix, x['var_chr'], x['var_from'], x['var_id'], f"{x['HostGenes']}:{x['HostGenes']}"), axis=1, result_type="expand")

NominalPass_df.to_csv("scratch/SpliceSiteEffects.txt.gz", sep='\t')
