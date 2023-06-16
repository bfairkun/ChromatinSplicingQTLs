#### Development rules written by Carlos. 
#### Some of these might be redundant with rule all, or with other rules
#### Others might be obsolete
        
rule SpliceSiteMutations:
    input:
        expand(
            "QTLs/QTLTools/{Phenotype}.{Prime}/NominalPass.txt.gz",
            Phenotype = ['chRNA.Splicing', 'polyA.Splicing', 'H3K4ME1', 'H3K4ME3', 'H3K27AC', 'H3K36ME3'],
            Prime = ['5PrimeSS', '3PrimeSS']
        ),
        "QTLs/QTLTools/polyA.Splicing.5PrimeSS.Subset_YRI/NominalPass.txt.gz",
        "QTLs/QTLTools/polyA.Splicing.3PrimeSS.Subset_YRI/NominalPass.txt.gz",
        

rule plotNonCodingRNA:
    input:
        expand("NonCodingRNA_annotation/bigwig/{IndID}.plus.bw", 
        IndID = ["NA18486", "NA18497", "NA18498", "NA18499", "NA19201", "NA19210"]),
        expand("NonCodingRNA_annotation/bigwig/{IndID}.minus.bw", 
        IndID = ["NA18486", "NA18497", "NA18498", "NA18499", "NA19201", "NA19210"]),
        expand("NonCodingRNA_annotation/deeptools/matrix/{IndID}.RNA{strictness}.matrix.gz", 
        IndID = ["NA18486", "NA19201", "NA19210"], strictness=["", ".strict1", ".strict2", ".strict3"]),
        expand("NonCodingRNA_annotation/deeptools/plots/{IndID}.RNA{strictness}.plot.png", 
        IndID = ["NA18486", "NA19201", "NA19210"], strictness=["", ".strict1", ".strict2", ".strict3"]),
        expand("NonCodingRNA_annotation/deeptools/matrix/{IndID}.histone{strictness}.matrix.gz", 
        IndID = ["NA19201", "NA19210"], strictness=["", ".strict1", ".strict2", ".strict3"]),
        expand("NonCodingRNA_annotation/deeptools/plots/{IndID}.histone{strictness}.plot.png",  
        IndID = ["NA19201", "NA19210"], strictness=["", ".strict1", ".strict2", ".strict3"]),
        expand("NonCodingRNA_annotation/deeptools/matrix/{IndID}.combined{strictness}.matrix.gz", 
        IndID = ["NA19201", "NA19210"], strictness=["", ".strict1", ".strict2", ".strict3"]),
        expand("NonCodingRNA_annotation/deeptools/plots/{IndID}.combined{strictness}.plot.png",  
        IndID = ["NA19201", "NA19210"], strictness=["", ".strict1", ".strict2", ".strict3"]),
        

  
rule ncRNA_annotation_caller:
    input:
        "NonCodingRNA_merged/annotation/ncRNA.bed.gz",
        'NonCodingRNA_merged/annotation/ncRNA.histone.tab.gz',
        'NonCodingRNA_merged/annotation/allGenes.histone.tab.gz',
        "NonCodingRNA_merged/Expression_ncRNA_and_reverse/chRNA.Expression_featureCounts/Counts.txt"
        
rule QTLTools_caller:
    input:
        expand(
          "QTLs/QTLTools/chRNA.Expression{RNA_type}/OnlyFirstReps.qqnorm.bed.gz",
          RNA_type = [".Splicing", "_ncRNA"],
        ),
        #expand("RPKM_tables/{Phenotype}.RPKM.bed.gz",
        #       Phenotype = ["chRNA", "polyA", "MetabolicLabelled.30min", "MetabolicLabelled.60min"]
        #),
        expand(
          "QTLs/QTLTools/{Phenotype}/{Pass}.txt.gz",
          Phenotype = ["chRNA.Splicing", "chRNA.Expression.Splicing", "chRNA.Expression_ncRNA",
                       "polyA.Splicing", "Expression.Splicing", "polyA.Expression_ncRNA", 
                       "polyA.Splicing.Subset_YRI", "Expression.Splicing.Subset_YRI", "polyA.Expression_ncRNA.Subset_YRI",
                       "MetabolicLabelled.30min_ncRNA", "MetabolicLabelled.60min_ncRNA",
                       "chRNA.Splicing.Order", "chRNA.IER", "polyA.IER", 
                       "polyA.IER.Subset_YRI", "chRNA.Slopes", "H3K36ME3", "H3K4ME3", "H3K4ME1","H3K27AC",
                       "H3K36ME3_ncRNA", "MetabolicLabelled.30min.IER", "MetabolicLabelled.60min.IER",                                "ProCap_uaRNA", "DNaseISensitivity"],
          Pass = ["PermutationPass.FDR_Added", "NominalPass"]
        ),


rule GatherCPM:
    input:
        expand(
          "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.qqnorm.bed.gz",
          Phenotype = ["chRNA.Expression.Splicing", "Expression.Splicing", 
                       "Expression.Splicing.Subset_YRI",  "MetabolicLabelled.30min", 
                       "MetabolicLabelled.60min", 'H3K36ME3'],
        ),
        expand(
        "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized_AtTSS.qqnorm.bed.gz",
        Phenotype = ["H3K27AC", "H3K4ME1", "H3K4ME3"]
        ),

rule HMM_develop:
    input:
        "NonCodingRNA/annotation/ncRNA.bed.gz",
        "NonCodingRNA/annotation/NonCodingRNA.bed.gz",
        "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
        "NonCodingRNA/annotation/allGenes.TSS.bed.gz",
        "NonCodingRNA/annotation/ncRNA.3PrimeEnd.bed.gz",
        "NonCodingRNA/annotation/tmp/uaRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/uaRNA.annotated.bed.gz",
        "NonCodingRNA/annotation/tmp/incRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/incRNA.annotated.bed.gz",
        "NonCodingRNA/annotation/ncRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/lncRNA.ncRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/srtRNA.bed.gz",
        'NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz',
        'NonCodingRNA/annotation/allGenes.annotation.tab.gz',
        
                
        
rule IntronCounts_develop:
    input:
        "SplicingAnalysis/GeneIntronCounts/chRNA.GeneIntronCounts.bed.gz",
        "SplicingAnalysis/GeneIntronCounts/polyA.GeneIntronCounts.bed.gz",
        "SplicingAnalysis/GeneIntronCounts/MetabolicLabelled.30min.GeneIntronCounts.bed.gz",
        "SplicingAnalysis/GeneIntronCounts/MetabolicLabelled.60min.GeneIntronCounts.bed.gz",
        

rule GatherDNaseBigWigs:
    input:
        expand(
            "bigwigs/DNaseISensitivity/{IndID}.merged.bw",
            zip,
            IndID=DNaseSamples_df["IndID"],
        ),


rule GatherProCapBigWigs:
    input:
        expand(
            "bigwigs/ProCap_stranded/{IndID}.{Rep}.plus.bw",
            zip,
            IndID=ProCapSamples_df["IndID"],
            Rep=ProCapSamples_df["RepNumber"],
        ),
        expand(
            "bigwigs/ProCap_stranded/{IndID}.{Rep}.minus.bw",
            zip,
            IndID=ProCapSamples_df["IndID"],
            Rep=ProCapSamples_df["RepNumber"],
        ),


rule MatrixForMetaplots:
    input:
        expand('NonCodingRNA/QTLs/metaplots/mat/{QTLs}.{var}.{phe}.mat.gz',
            QTLs = ['diQTLs', 'apQTLs'],
            var = ['ncQTL', 'eQTL'],
            phe = ['uaRNA', 'uaGene']),
        expand('NonCodingRNA/QTLs/metaplots/mat/{QTLs}.{var}.{phe}_metagene.mat.gz',
            QTLs = ['diQTLs', 'apQTLs'],
            var = ['ncQTL', 'eQTL'],
            phe = ['uaGene'])
            
            
rule QTLTools_500kb:
    input:
        expand(
          "QTLs/QTLTools/{Phenotype}/{Pass}.txt.gz",
          Phenotype = ["H3K36ME3", "H3K4ME3", "H3K4ME1","H3K27AC"],
          Pass = ["PermutationPass500kb.FDR_Added", "NominalPass500kb"]
        ),
        
rule QTLTools_250kb:
    input:
        expand(
          "QTLs/QTLTools/{Phenotype}/{Pass}.txt.gz",
          Phenotype = ["H3K36ME3", "H3K4ME3", "H3K4ME1","H3K27AC"],
          Pass = ["PermutationPass250kb.FDR_Added"]
        ),

rule GatherLongReads:
    input:
        expand(
          "LongReads/Junctions/{sample}.junc.filtered.bed.gz",
          sample = long_read_samples
        ),
          
          
rule GatherONTFastq:
    input:
        expand("Fastq/ONT/{Phenotype}/{IndID}/R1.fastq.gz", zip, Phenotype=long_read_samples_df['Phenotype'], IndID=long_read_samples_df['IndID'])
        
rule GatherLongReadJunc:
    input:
        expand("LongReads/Junctions/{sample}.annotated.junc.gz", sample = all_long_reads),
        #'LongReads/Analysis/IsoSeq.nmd.tab.gz',
        #'LongReads/Analysis/IsoSeq.stable.tab.gz',
        #'LongReads/Analysis/IsoSeq.nmd_avg.tab.gz',
        #'LongReads/Analysis/IsoSeq.stable_avg.tab.gz',
        #'LongReads/Analysis/NMD_KD.nmd.tab.gz',
        #'LongReads/Analysis/NMD_KD.stable.tab.gz',
        #'LongReads/Analysis/NMD_KD.nmd_avg.tab.gz',
        #'LongReads/Analysis/NMD_KD.stable_avg.tab.gz',
        #'LongReads/Analysis/Churchman.nmd.tab.gz',
        #'LongReads/Analysis/Churchman.stable.tab.gz',
        #'LongReads/Analysis/Churchman.nmd_avg.tab.gz',
        #'LongReads/Analysis/Churchman.stable_avg.tab.gz',
        #'LongReads/Analysis/IsoSeq.ByQuartile.tab.gz',
        #'LongReads/Analysis/NMD_KD.ByQuartile.tab.gz',
        #'LongReads/Analysis/Churchman.ByQuartile.tab.gz'
        


rule NMDJunctions_test:
    input:
        expand(
            "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.{ext}",
            Junction = ["FivePrime", "ThreePrime", "ThreePrime.ForMaxEntScan"], 
            ext = ['tab', 'fa']
        ),
        expand(
            "SplicingAnalysis/NMDJunctions/tables/leafcutter_perind.counts_only.{Phenotype}.bed.gz",
            Phenotype = ["Expression.Splicing", "chRNA.Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min"]),
        expand("SplicingAnalysis/NMDJunctions/Annotation/{Junction}.scored.tab.gz",
        Junction = ["FivePrime", "ThreePrime"]),
        expand("SplicingAnalysis/NMDJunctions/Plots/{Junction}.PWM.png",
        Junction = ["FivePrime", "ThreePrime"]),
        expand("SplicingAnalysis/NMDJunctions/Annotation/{Junction}.MaxEntScan_score.tab.gz",
        Junction = ["FivePrime", "ThreePrime.ForMaxEntScan"]),
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored.bed.gz",
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_junctions.bed.gz",
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_clusters.bed.gz",
     
rule MakeEURPhenotypes:
    input:
        expand("QTLs/QTLTools/{Phenotype}.Subset_EUR/OnlyFirstReps.qqnorm.bed.gz",
            Phenotype = ["Expression.Splicing", "polyA.Splicing"]),
        expand("QTLs/QTLTools/{Phenotype}.Subset_EUR/{Pass}.txt.gz",
            Phenotype = ["Expression.Splicing", "polyA.Splicing"],
            Pass = ["NominalPass", "PermutationPass.FDR_Added"]),
        
rule CollectBasicMapping:
    input:
        expand("featureCountsBasicGtf/{Phenotype}/Counts.txt", 
        Phenotype=["chRNA.Expression.Splicing","Expression.Splicing","MetabolicLabelled.60min", "MetabolicLabelled.30min"]),
        expand("QTLs/QTLTools/{Phenotype}/OnlyFirstRepsBasicUnstandardized.qqnorm.bed.gz",
        Phenotype=["chRNA.Expression.Splicing","Expression.Splicing","MetabolicLabelled.60min", "MetabolicLabelled.30min"]),
        "ReadLengthMapExperimentResults/tables/AllRNASeq.Normalized.RPKM.bed.gz",
        "QTLs/QTLTools/H3K4ME3/CountsPeaksAtTSS.bed.gz",
        "QTLs/QTLTools/H3K27AC/CountsPeaksAtTSS.bed.gz",
        
        
rule CollectFinemapping:
    input:
        expand("FineMapping/susie_runs_{Subset}/susie_output.tab.gz",
               Subset = ["Geuvadis", "EUR", "YRI"]),
        expand("FineMapping/susie_runs_{Subset}/susie_output_long.tab.gz",
               Subset = ["Geuvadis", "EUR", "YRI"]),
        #"FineMapping/Genotypes/1KG_GRCh38/EUR.txt.bgz.tbi",
        #"FineMapping/susie_runs/susie_output.tab"
          
rule collect_metaplot:
    input:
        #expand("Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{metaplot}.png",
        #    Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min",
        #    "MetabolicLabelled.60min", "H3K36ME3"], IndID = ['NA19137', 'NA19152', 'NA19153'],
        #    metaplot=["CDS", "gene"]),
        expand("Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{quartiles}.{stat}.png",
            Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min",
            "MetabolicLabelled.60min", "H3K36ME3", "H3K27AC", "H3K4ME3", "H3K4ME1"], 
            IndID = ['NA19137', 'NA19152', 'NA19153'], stat = ['mean', 'median'],
            quartiles = [#"CDS_length", "gene_length", 
            "expression"]),
        #expand("Metaplots/AssayProfiles/Plots/chRNA.Expression.Splicing.{IndID}.{metaplot}.{strand}.png",
        #    IndID = ['NA18486', 'NA19137', 'NA19152', 'NA19153'],
        #    metaplot=["CDS", "gene"], strand = ["plus", "minus"]),
        expand("Metaplots/AssayProfiles/Plots/chRNA.Expression.Splicing.{IndID}.{quartiles}.{strand}.{stat}.png",
            IndID = ['NA18486', 'NA19137', 'NA19152', 'NA19153'], stat = ['mean', 'median'],
            quartiles = [#"CDS_length", "gene_length", 
            "expression"], strand = ["plus", "minus"]),
        "QTLs/QTLTools/H3K36ME3/OnlyFirstRepsUnstandardized.Last3K.qqnorm.bed.gz",
        "QTLs/QTLTools/H3K36ME3/OnlyFirstRepsUnstandardized.Encode.qqnorm.bed.gz",
        expand("Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{quartiles}.{stat}.png",
            Phenotype = ["H3K36ME3"], IndID = ['NA18486'], stat = ['mean', 'median'],
            quartiles = ["expression"]),#, "CDS_length", "gene_length"]),
        #expand("Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{metaplot}.png",
        #    Phenotype = ["H3K36ME3"], IndID = ['NA18486'],
        #    metaplot=["CDS", "gene"]),
            
rule CollectMetaplotIntrons:
    input:
        expand("Metaplots/AssayProfiles/Plots/Introns.{IndID}.{strand}.png", 
        IndID = ['NA18486', 'NA19137', 'NA19152', 'NA19153'], strand = ['plus', 'minus']),
        expand("Metaplots/AssayProfiles/Plots/Introns.{IndID}.{strand}.long_introns.png", 
        IndID = ['NA18486', 'NA19137', 'NA19152', 'NA19153'], strand = ['plus', 'minus'])
            
rule H3K36ME3_corr:
    input:
        "QTLs/QTLTools/H3K36ME3/OnlyFirstRepsUnstandardized.Encode.qqnorm.bed.gz",
        "QTLs/QTLTools/H3K36ME3/OnlyFirstRepsUnstandardized.qqnorm.bed.gz",
        "QTLs/QTLTools/H3K36ME3/OnlyFirstRepsUnstandardized.Last3K.qqnorm.bed.gz"
        
rule CollectTES_TSS:
    input:
        expand("LongReads/ClosestAnnotatedTermini/{GTFType}_{TES_or_TSS}/{sample}.bed.gz",
               GTFType = ["GTFTools_BasicAnnotations", "GTFTools"], TES_or_TSS = ["TSS", "TES"],
               sample = all_long_reads)

