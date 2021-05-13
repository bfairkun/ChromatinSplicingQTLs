# Chromatin, splicing and genetic variation project

Planning to integrate published and new data to explore lots of ideas related to chromatin patterns, co-txnal splicing.

### Questions of interest:
- Basic characterization of stuff we can do without considering variation between individuals
- what sequence features can regulate H3K36me3, RNAPII-phosphorylation. Much literature about how these things are established over gene bodies, H3K36me <--> RNAPIISer2P interactions, Ser2 phosphorylation and H3K36me3 accumulates as transcription proceeds to 3' end... But not much literature on the sequence features that influence these patterns.
- How does splicing influence chromatin and transcription. Quantifying H3K4me3 and TSSs could be of interest [E-mats mechanism](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5098600/). Can we see this effect. Characterize distance constraints of mechanism
- H3K9me3 usual model with regards to splicing is that it slows polymerase, promotes exon inclusion. So do colocalized sQTL/H3K9me3QTL/polymeraseQTL generally have concordant effect directions with that. Same with H3K36me3
- Can we anchor some splicing/transcription/chromatin causal graphs for colocolized QTLs based on well-mapped QTL positions: Like if a QTL is in a TF binding site in an enhancer, versus if it is in a splice site. What are directions of effects (Like [McVicker et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3947669/))
- Measuring variants that effect splicing/transcription speed. Do they associate with PolII phosphorylation?
- Assess contribution of any new molecular phenotypes to complex traits
- Are there any GWAS variants that we can see as sQTL/eQTL in chRNA but not polyARNA. That is interesting in itself I think. Presumably GWAS variants work closer to cytoplasm, so maybe these still have consistent effects in polyA fraction that are just smaller or harder to estimate due to NMD or other rapid decay mechanisms.

### Experiments overview/plans:
- ~50 fractionations done so far, in batches of 20-30 each for cell growth and fractionation, then a batch of 48 for chRNA-seq library prep.
- Currently growing cells for another fractionation batch. Planning one more batch after that (finish in a month-ish).
- Library prep for chRNA-seq + Cut&Tag another couple weeks after that. So long as the preliminary Cut&Tag data look good, I think we should have all the data sequenced sometime in July.

### Data processing
- see [code/config/samples.tsv](https://github.com/bfairkun/ChromatinSplicingQTLs/blob/master/code/config/samples.tsv) for full list of published and in-house samples that I have already started processing.

#### Data sources:
Will have data for all the following assays for YRI panel (and geuvadis has europeans too)
##### Stuff I have already started processing, or writing code to process
- H3K4me3, H3K4me1, H3K27Ac ChIP-seq: [Grubert et al](https://pubmed.ncbi.nlm.nih.gov/26300125/)
- RNA-seq: geuvadis
- chRNA-seq: in-house
- H3K36me3, RNAPII-Ser5P, RNAPII-Ser2P, possibly H3K79me2 & H3K9me3 Cut&Tag: in house (still need to check data quality from tests, but will probably prep libraries for all of those marks for YRI panel)
- CTCF ChIP-seq: [Ding et al](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004798)
##### Published Stuff I haven't done anything for yet, but I am aware of that could be relevant
- Could also consider pro-cap/TSS initiation: [Kristjánsdóttir et al](https://www.nature.com/articles/s41467-020-19829-z)
- 4sU, RNA methylation, DNA methylation, DNAseI

#### basic NGS processing:
- adapter clipping (fastp)
- alignment w/ STAR (RNA-seq) or hisat2 (ChromatinProfiling), used WASP filtering
- QC: multiqc, qualimap (check that chRNA maps intron mapping rate)

#### Stuff to look at without focusing on genetic variation
- chRNA characterization: meta-intron coverage. Enrichment for novel junctions, NMD-specific junctions
- Obtaining genome-wide splicing timing estimate from the intron slope thing (something I showed you a preliminary analysis for) now that I think of it is probably too complicated because of unknown amount of fragmentation. More fragmentation before library prep will lead to a flatter slope for non-biological reasons, biasing estimates in ways that are hard to reasonably correct for. This shouldn't really matter for QTL calling, since we are just looking for relative differences by genotype.
- eRNA representation in chRNA-seq
- co-transcriptional splicing is more efficient for constitutive introns than ones flanking cassette exons (previously I showed this using [3SS ratio method](https://www-sciencedirect-com.proxy.uchicago.edu/science/article/pii/S1046202315001711) (Normalized to polyA-RNA-seq) for quantifying co-txnal splicing efficiency). Based on the much research on splicing efficiency correlating with highly active transcription, I would expect to see some correlation with gene expression too, which I haven't yet looked for.

#### Phenotype quantification

- ChromatinProfiling peak calling to identify phenotype features for testing with macs2... Testing peak calls with all individuals jointly, or with each individual and then merging.
- For RNA-seq, top 14000 protein coding genes (as measured by median CPM in geuvadis) counted with featureCounts and standardized and inverse-normalized to association testing.
- For ChIP-seq, same process with the macs2 peak calls.
- For H3K36me, (don't have data yet) I may test out a couple ideas for determining test features. Like just using the same 14000 expressed gene regions, and then breaking those regions into 1000bp windows
- Splicing junctions quantified with leafcutter pipeline.
- Intron slope thing that I showed before. I wonder how much of those slope affects are more related to transcription initiation than elongation rate.

#### association testing
- QTL testing done with QTLtools, grouping splicing phenotypes by cluster to get Cluster-level tests (if I turn out breaking broad H3K36me3 phenotypes into windows, will group these too).
- Will include 0, 1, 2 or 3 genotype PCs as covariates, after manually inspecting PC plots for stratification. Will include expression PCs too: Planning to include the number of PCs for which the variance explained by a PC exceeds the variance explained from PCA on a permuted phenotype matrix (each row randomly shuffled). This is what [Rasqual paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5098600/) did, and from my testing this seems to give a reasonable number of PCs.
- For intron retention I will try simple approach: What makes sense to me is just testing constitutive introns (that don't overlap any annotated exons), then summing unspliced read counts over the intron, to quantify, and then filtering for well-expressed introns that are also efficiently spliced in standard RNA-seq, and standardizing and inverse-normalizing phenotypes like above before testing.
- For chromatin phenotypes I could consider doing the rasqual pipeline which I have working, I'll have to see how much more power it gets.


#### colocalization
- testing a lot of the questions will involve looking for colocalizing molQTLs. I'm not yet sure of how best to do this with all the different assays we have. could do this formally using coloc software in pairwise combinations? I know a multi-trait version of coloc exists as well. thoughts?

