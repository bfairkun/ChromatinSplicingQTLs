# liftOverBedpe

This is a simple python wrapper for the tool liftOVer that allows users to liftOver BEDPE files from one genome build to another.

## Getting Started

Just download and run.

### Prerequisites

liftOverBedpe requires [liftOver](https://genome.sph.umich.edu/wiki/LiftOver), a liftOver chain file, and python

### Installing

None required

### Example

```
# download loopcalls
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71831/suppl/GSE71831_Patski_HiCCUPS_looplist.txt.gz

# unzip loopcalls
gunzip GSE71831_Patski_HiCCUPS_looplist.txt.gz 

# format to BEDPE
cat GSE71831_Patski_HiCCUPS_looplist.txt | tail -n +2 | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,"loop_"NR,$16,"+","-"}' >  GSE71831_Patski_HiCCUPS_looplist.mm10.bedpe

# download liftOver chain
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz

# check out the help on the liftOverBedpe.py script (The exact path to the script might change depending on where it is on your computer)
python ../Tools/liftOverBedpe/liftOverBedpe.py -h

# convert to mm9 (The exact path to files depend on where files are on your computer)
python ../Tools/liftOverBedpe/liftOverBedpe.py --lift ~/Tools/liftOver --chain mm10ToMm9.over.chain.gz --i GSE71831_Patski_HiCCUPS_looplist.mm10.bedpe --o GSE71831_Patski_HiCCUPS_looplist.mm9.bedpe  --h F
```
