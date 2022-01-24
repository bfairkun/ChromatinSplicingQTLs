import pandas as pd

intron_bed = pd.read_csv('Misc/GencodeHg38_all_introns.corrected.uniq.bed', sep='\t',
                   names = ['chrom', 'start', 'end', 'names', 'score', 'strand'])

intron_SAF = pd.DataFrame()
intron_SAF['GeneID'] = list(intron_bed.names)
intron_SAF['Chr'] = list(intron_bed.chrom)
intron_SAF['Start'] = list(intron_bed.start+1)
intron_SAF['End'] = list(intron_bed.end)
intron_SAF['Strand'] = list(intron_bed.strand)

intron_SAF.to_csv('../data/introns.saf', index=False, header=True, sep='\t')
