import numpy as np
import pandas as pd

def longtable_by_annotation(long_table):
    counts_by_annotation = long_table.groupby(['Dataset', 'IndID', 'SuperAnnotation']).Count.sum().reset_index()
    counts_by_semiannotation = long_table.groupby(
        ['Dataset', 'IndID', 'SemiSupergroupAnnotations']
    ).Count.sum().reset_index()
    return counts_by_annotation, counts_by_semiannotation

def get_junction_types(counts_by_annotation, counts_by_semiannotation):
    
    counts = counts_by_annotation.groupby(['Dataset', 'IndID']).Count.sum().reset_index()
    counts.columns = ['Dataset', 'IndID', 'total_counts']

    assay = []
    for i in counts.Dataset:
        if i == 'chRNA.Expression.Splicing':
            assay.append('chRNA')
        elif i == 'Expression.Splicing':
            assay.append('polyA')
        elif i == 'MetabolicLabelled.30min':
            assay.append('4sU 30min')
        elif i == 'MetabolicLabelled.60min':
            assay.append('4sU 60min')

    counts['assay'] = assay

    Protein_coding_counts = np.array(counts_by_annotation.loc[counts_by_annotation.SuperAnnotation.isin(
        ['AnnotatedJunc_UnproductiveCodingGene',
         'UnannotatedJunc_UnproductiveCodingGene',
         'AnnotatedJunc_ProductiveCodingGene',
         'UnannotatedJunc_ProductiveCodingGene'])].groupby(['Dataset', 'IndID']).Count.sum())

    Unproductive_counts = np.array(counts_by_annotation.loc[counts_by_annotation.SuperAnnotation.isin(
        ['AnnotatedJunc_UnproductiveCodingGene',
         'UnannotatedJunc_UnproductiveCodingGene'])].groupby(['Dataset', 'IndID']).Count.sum())

    Productive_counts = np.array(counts_by_annotation.loc[counts_by_annotation.SuperAnnotation.isin(
        ['AnnotatedJunc_ProductiveCodingGene',
         'UnannotatedJunc_ProductiveCodingGene'])].groupby(['Dataset', 'IndID']).Count.sum())
    
    Annotated_NMD_counts = np.array(counts_by_semiannotation.loc[counts_by_semiannotation.SemiSupergroupAnnotations.isin(
        ['uniquely nonsense_mediated_decay tag'])].groupby(['Dataset', 'IndID']).Count.sum())
    
    unprod_annotations = ['uniquely processed_transcript tag',
       'overlaps processed transcript intron',
       'uniquely retained_intron tag', 'predicted_NMD pstopcodon',
       'predicted_NMD far3p', 'uniquely nonsense_mediated_decay tag',
       'predicted_NMD UTRjunction', 'predicted_NMD YN',
       'predicted_NMD far5p', 'predicted_NMD reason2',
       'overlaps retained_intron tag', 'stable.NY',
       'overlaps nonsense_mediated_decay intron',
       'uniquely non_stop_decay tag', 'predicted_NMD reason1',
       'predicted_NMD NN']

    predicted_annot = [x for x in unprod_annotations if ('predicted_NMD' in x)]

    try:
        Unannotated_NMD_counts = np.array(counts_by_semiannotation.loc[counts_by_semiannotation.SemiSupergroupAnnotations.isin(
            predicted_annot)].groupby(['Dataset', 'IndID']).Count.sum())

        counts['Unannotated_NMD'] = 100*(Unannotated_NMD_counts/Protein_coding_counts)
    except:
        print('No unannotated NMD')

    Annotated_pt_counts = np.array(counts_by_semiannotation.loc[counts_by_semiannotation.SemiSupergroupAnnotations.isin(
        ['uniquely processed_transcript tag'])].groupby(['Dataset', 'IndID']).Count.sum())
    
    counts['Unproductive_juncs'] = 100*(Unproductive_counts/Protein_coding_counts)
    counts['Productive_juncs'] = 100*(Productive_counts/Protein_coding_counts)

    counts['Annotated_NMD'] = 100*(Annotated_NMD_counts/Protein_coding_counts)
    
    counts['Processed_transcripts'] = 100*(Annotated_pt_counts/Protein_coding_counts)
    
    other_unproductive_annot = [x for x in unprod_annotations if (x not in [
        'uniquely processed_transcript tag', 'uniquely nonsense_mediated_decay tag'] + predicted_annot)]
    
    other_unproductive = np.array(counts_by_semiannotation.loc[counts_by_semiannotation.SemiSupergroupAnnotations.isin(
            other_unproductive_annot)].groupby(['Dataset', 'IndID']).Count.sum())
    
    counts['other_unproductive'] = 100*(other_unproductive/Protein_coding_counts)
    
    

    order = ['chRNA', '4sU 30min', '4sU 60min', 'polyA']
    
    return counts, order

def longtable_to_boxplot(long_table):
    
    counts_by_annotation, counts_by_semiannotation = longtable_by_annotation(long_table)
    counts, order = get_junction_types(counts_by_annotation, counts_by_semiannotation)
    
    return counts, order