#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PlotQTLs
# @created     : Wednesday Apr 20, 2022 15:38:03 CDT
#
# @description : 
######################################################################

import sys
import pandas as pd
sys.path.insert(1, 'scripts/GenometracksByGenotype')
import AggregateBigwigsForPlotting as qtlplot
import subprocess
import numpy as np
import time

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "PlotQTLs/Temp/", "scratch/testplot.", "1" ,"500"]

temp_prefix_no_chunk  = sys.argv[1]
final_output_prefix = sys.argv[2]
chunk_n, TotalChunks_N = [int(i) for i in sys.argv[3:5]]

temp_prefix = f'{temp_prefix_no_chunk}chunk{chunk_n}.'
temp_bed_fn = f'{temp_prefix}ColocalizedFeats.bed'

dat = pd.read_csv("../output/hyprcoloc_results/ForColoc/MolColocStandard/hyprcoloc.results.OnlyColocalized.Stats.txt.gz", sep='\t')
dat['score'] = dat['Locus'] + '_' + dat['snp']


all_feats = pd.read_csv("PlotQTLs/ColocTestFeatures.bed.gz", sep='\t', header=0, names=["Chr", "start", "stop", "name", "score", "strand", "thickStart", "thickEnd", "color"])
coloc_feats = pd.read_csv("PlotQTLs/ColocFeatures.bed.gz", sep='\t', header=0, names=["Chr", "start", "stop", "name", "score", "strand", "thickStart", "thickEnd", "color"])

ChrSizes = pd.read_csv("ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai", usecols=[0,1], sep='\t', names=["Chr", "size"])

joined = dat.merge(coloc_feats, how='inner', on='score').merge(ChrSizes, how='inner', on='Chr')
joined['max'] = joined.groupby(['score'])['stop'].transform('max')
joined['min'] = joined.groupby(['score'])['start'].transform('min')
joined['range'] = joined['max'] - joined['min']
joined['left'] = joined.apply (lambda row: max([row['min'] - row['range']*0.1, 1]), axis=1).astype(int)
joined['right'] = joined.apply (lambda row: min([row['max'] + row['range']*0.1, row['size']]), axis=1).astype(int)
joined['DummyScore'] = 0

ToPlot = joined[['snp', 'Chr', 'left', 'right', 'score', 'DummyScore']].drop_duplicates()
ToPlot['snpPos'] = ToPlot.apply (lambda row: row['Chr'] + ':' + row['snp'].split(':')[1], axis=1)
# ToPlotSubset = ToPlot.sample(10)
ToPlotSubset = ToPlot.iloc[np.array_split(range(len(ToPlot.index)), TotalChunks_N)[chunk_n-1]]

print(ToPlotSubset)

for i, row in ToPlotSubset.iterrows():
    try:
        joined.loc[joined['score'] == row['score']][['Chr', 'start', 'stop', 'name', 'DummyScore' ,'strand', 'thickStart', 'thickEnd', 'color']].to_csv(temp_bed_fn, sep='\t', index=False, header=False)
        MyArgs = f"--Region {row['Chr']}:{row['left']}-{row['right']} --SnpPos {row['snpPos']} --VCF Genotypes/1KG_GRCh38/{row['Chr'][3:]}.vcf.gz --BigwigListType KeyFile --GroupSettingsFile PlotQTLs/bwList.Groups.tsv --BigwigList PlotQTLs/bwList.tsv --OutputPrefix {temp_prefix} --FilterJuncsByBed {temp_bed_fn} --TracksTemplate scripts/GenometracksByGenotype/tracks_templates/GeneralPurposeColoredByGenotype.ini"
        parsed_args = qtlplot.parse_args(filter(None, MyArgs.split(' ')))
        qtlplot.main(**vars(parsed_args))
        with open(f"{temp_prefix}tracks.ini", 'a') as tracks_ini_fh:
            _ = tracks_ini_fh.write(
                    f'''
                    [test_feats]
                    title = test_feats
                    height = 2
                    # gene_rows = 10
                    fontsize = 4
                    file_type = bed
                    file = PlotQTLs/ColocTestFeatures.bed.gz
                    color = bed_rgb
                    color_backbone = bed_rgb
                    all_labels_inside = true
                    labels_in_margin = true
                    [line]
                    file_type = hlines
                    y_values = 0.5
                    show_data_range = false
                    [colocalized_feats]
                    title = colocalized_feats
                    height = 2
                    # gene_rows = 10
                    fontsize = 4
                    file_type = bed
                    file = {temp_bed_fn}
                    color = bed_rgb
                    color_backbone = bed_rgb
                    all_labels_inside = true
                    labels_in_margin = true
                    [x-axis]
                    [spacer]
                    height = 0.5
                    [genes]
                    file = scripts/GenometracksByGenotype/PremadeTracks/gencode.v26.FromGTEx.genes.bed12.gz
                    height = 1
                    color = darkblue
                    labels = true
                    fontsize = 10
                    style = UCSC
                    file_type = bed
                    all_labels_inside = true
                    labels_in_margin = true
                    '''
                    )
        subprocess.run(f"pyGenomeTracks --region {row['Chr']}:{row['left']}-{row['right']} --tracks {temp_prefix}tracks.ini -out {final_output_prefix}{row['score'].replace(':', '_')}.pdf -t {row['score'].replace(':', '_')} --trackLabelFraction 0.08".split(' '))
    except:
        "Problem! Skipping..."


