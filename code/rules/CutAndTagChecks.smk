DownloadEncode_FC_Over_Input_ChIPSeq_bw_dict = {
        "POL2S5" : "https://www.encodeproject.org/files/ENCFF559YIO/@@download/ENCFF559YIO.bigWig",
        "POL2S2" : "https://www.encodeproject.org/files/ENCFF353DID/@@download/ENCFF353DID.bigWig",
        "H3K79ME2" : "https://www.encodeproject.org/files/ENCFF931USZ/@@download/ENCFF931USZ.bigWig",
        "H3K36ME3" : "https://www.encodeproject.org/files/ENCFF380LZI/@@download/ENCFF380LZI.bigWig",
        "H3K9ME3" : "https://www.encodeproject.org/files/ENCFF683HCZ/@@download/ENCFF683HCZ.bigWig"
        }

localrules: GatherEncodeBws, DownloadEncode_FC_Over_Input_ChIPSeq_bw

rule GatherEncodeBws:
    input:
        expand("ENCODE/ChIPSeq/FC_Over_input/{Antibody}.bw", Antibody=DownloadEncode_FC_Over_Input_ChIPSeq_bw_dict.keys())

rule DownloadEncode_FC_Over_Input_ChIPSeq_bw:
    output:
        "ENCODE/ChIPSeq/FC_Over_input/{Antibody}.bw"
    params:
        lambda wildcards: DownloadEncode_FC_Over_Input_ChIPSeq_bw_dict[wildcards.Antibody]
    shell:
        "wget -O {output} {params} "

use rule MakeBigwigs as MakeBigwigs_Unfiltered with:
    input:
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        bam = GetUnfilteredBamForBigwig,
        bai = GetUnfilteredBaiForBigwig
    log:
        "logs/MakeBigwigs_Unfiltered/{Phenotype}/{IndID}.{Rep}.bw"
    output:
        bw = "bigwigs_FromNonWASPFilteredReads/{Phenotype}/{IndID}.{Rep}.bw"
