# RNA-hyper-editing

Our README walks users through how to utilize all the tools in this repository to detect RNA hyper-editing in RNA-sequencing data sets, expanding on the methodology outlined in [Cuddleston et al. 2021](https://www.biorxiv.org/content/10.1101/2021.06.10.447947v1).

This pipeline was adapted from [Porath et al. 2014](https://www.nature.com/articles/ncomms5726). The original scripts are publicly available on their [github](https://github.com/hagitpt/Hyper-editing). The changes we have made are outlined below and also clearly noted in comments throughout.

## Directory Structure

```
|-experiment_name
  |-GenomeIndex
  |-GenomeIndexTrans
  |-H276_GABA
    |-HE_scripts
```
Inside of the master directory for your experiment, you create an output directory for each sample. This is specified in the "MODIFY PER SAMPLE" section of config.sh (line 11). In each sample output directory, recursively place the HE_scripts directory, which contains all of the scripts needed for hyper-editing detection.

## Indexing and Transforming the Reference Genome
The script TransformIndexBWA_genome.sh requires a genome reference file provided as input with ".fa" extension, however the file extension must be left out when executing from the command line. An example is provided below for the genome file named GRCh38.primary_assembly.genome.fa as input:
```unix
TransformIndexBWA_genome.sh GRCh38.primary_assembly.genome
```
The path to the original geonome and transformed genome indices should be provided in the "HARD CODED" section of the config.sh (lines 14-16).

## Config File
#### HARD CODED
Provide the full paths to the GenomeIndex and GenomeIndexTrans directories. Specify whether input is paired-end (1) or single-end (0). Here you can also modify the PHRED score offset, maximum gap size between pairs (if data is paired-end), specify to have the pipeline skip the BWA alignment step, and specify other options. Defaults are provided in this example.

#### MODIFY PER SAMPLE
On line 9 we provide a command to generate the file list on the fly for each sample individually, rather than having one large list for every sample in the experiment. However the file list is generated, provide the full path for it and also specify the output directory for each sample.

#### DETAILS FROM THE CONFIG FILE
Provide full paths to BWA and Samtools software.

## runHyper File
#### HARD CODED
In the sections "scripts" and "files and directories", provide the full file path to the HE_scripts housed individually within each sample output directory.

A few places throughout the runHyper script we found it to be helpful to hardcode the path, either to components of the hyper-editing pipeline (lines 48-49) or to the unmapped bam files (lines 79, 82, and 85). These places are also notated with comments throughout the runHyper script itself.

## Additional Processing Steps
#### Step One: Create UE.bed and ES.bed Files
Generate two files in bed format from the output of the pipeline, one with the coordinates of the A2G hyper-edited clusters (UE.bed) and one with the coordinates of the A2G sites (ES.bed). Since the example we have been providing used paired-end data, the pipeline will generate two separate output files which must be concatenated together.

```unix
PATH1="/full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-1.UE.bed_files"
PATH2="full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-2.UE.bed_files"
cat ${PATH1}/A2G.bed ${PATH2}/A2G.bed > H276_GABA.UE.bed
```

```unix
PATH1="/full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-1.ES.bed_files"
PATH2="full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-2.ES.bed_files"
cat ${PATH1}/A2G.bed ${PATH2}/A2G.bed > H276_GABA.ES.bed
```

#### Step Two: Merging Overlapping Hyper-Editing Clusters and Counting Number of Editing Sites Per Cluster
First we merge the hyper-edited clusters, so that any objects with overlapping coordiantes would be reported as only one hyper-editing event and not two separate clusters. Then we use the bed file of the A2G editing events within the hyper-editing clusters to count the number of editing sites per cluster, which is used in the Step Three of filtering. The following lines of code require [Bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html).

```unix
sort -k1,1 -k2,2n H276_GABA.UE.bed > H276_GABA.UE.bed.sorted

bedtools merge -i H276_GABA.UE.bed.sorted > H276_GABA.UE.bed.sorted.merged

bedtools intersect -a H276_GABA.UE.bed.sorted.merged -b H276_GABA.ES.bed -wa -c > H276_GABA.UE.bed.counted.temp
```

###  Step Three: Stretching Hyper-Edited Clusters
On a sample level, the average length (in nucleotides) of each hyper-editing cluster is calculated, in addition to the average number of A2G editing sites within these hyper-editing clusters. This number is then subtracted from the "start" coordinate and added to the "end" coordinate, effectively lengthening the cluster by 2x the average distance between editing sites, in order to account for the possibility that low-sequencing coverage at the edges of the defined cluster borders prevents inclusion of editing events resonably nearby.

```R
UE <- read.delim("H276_GABA.UE.bed.counted.temp", header = FALSE, sep = "\t")
totalclusters <- nrow(UE)

len <- UE$V3-UE$V2
avg.cluster.len <- mean(len)
ES <- UE$V4
avg.ES.per.cluster <- mean(ES)
ESdis <- avg.cluster.len/avg.ES.per.cluster
ESdis <- round(ESdis)
X <- type.convert(ESdis)
chr <- as.vector(UE$V1, mode = "character")
start <- as.vector(UE$V2, mode = "numeric")
end <- as.vector(UE$V3, mode = "numeric")

new.start <- UE$V2 - X
new.end <- UE$V3 + X

UEstretched <- cbind(chr, new.start, new.end)
UEstretched <- as.data.frame(UEstretched)
write.table(UEstretched, file = "H276_GABA.UE.stretched", sep = "\t", append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)

```

#### Step Four: Merging Overlapping Clusters and Counting Number of Editing Sites Per Cluster One Final Time
After adjusting the boundaries of the cluter, or "stretching" as we refer to it, merging overlapping clusters and counting the number of A2G editing events in each hyper-editing cluster is performed one final time.

```unix
bedtools merge -i H276_GABA.UE.stretched > H276_GABA.UE.bed.final

bedtools intersect -a H276_GABA.UE.bed.final -b H276_GABA.ES.bed -wa -c > H276_GABA.UE.bed.counted.final

```
