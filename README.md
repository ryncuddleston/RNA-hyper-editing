# Hyperediting Detection Pipeline

Our README walks users through how to utilize all the tools in this repository to detect RNA hyper-editing in RNA-seq datasets, expanding on the methodology outlined in [Cuddleston & Li et al., 2022](https://www.nature.com/articles/s41467-022-30531-0)

This pipeline was adapted from [Porath et al. 2014](https://www.nature.com/articles/ncomms5726). The scripts are publicly available on their [github](https://github.com/hagitpt/Hyper-editing). The changes we've made are outlined below and clearly noted in comments throughout.

## Directory Structure
Inside of the master directory for your experiment, you create an output directory for each sample. This is specified in the "MODIFY PER SAMPLE" section of config.sh (line 11). In each sample output directory, recursively place the HE_scripts directory, which contains all of the scripts needed for hyper-editing detection.

```
|-experiment_name
  |-GenomeIndex
  |-GenomeIndexTrans
  |-H276_GABA
    |-HE_scripts
```

Under the branch experiment_name branch we provide the unmapped fastq files, output following alignment with [STAR](https://github.com/alexdobin/STAR) for users to work with as example data for all the code below. 

## Indexing and Transforming the Reference Genome
The script TransformIndexBWA_genome.sh requires a genome reference file provided as input with ".fa" extension, however the file extension must be left out when executing from the command line. An example is provided below for the genome file named GRCh38.primary_assembly.genome.fa as input:

```
TransformIndexBWA_genome.sh GRCh38.primary_assembly.genome
```

The path to the original genome and transformed genome indices needs to be provided into the "HARD CODED" section of config.sh (lines 14-16).

## Config File
#### HARD CODED
Provide the full paths to the GenomeIndex and GenomeIndexTrans directories. Specify whether input is paired-end (1) or single-end (0). Here you can also modify the PHRED score offset, maximum gap size between pairs (if data is paired-end), specify to have the pipeline skip the BWA alignment step, and specify other options. Default parameters are used in this example.

#### MODIFY PER SAMPLE
On line 9 we provide a command to generate the file list on the fly for each sample individually, rather than having one large list for every sample in the experiment. This facilitates parallel execution of analyses. However the file list is generated, provide the full path for it and also specify the output directory for each sample.

#### DETAILS FROM THE CONFIG FILE
Provide full paths to BWA and Samtools software.

## runHyper File
#### HARD CODED
In the sections "scripts" and "files and directories", provide the full file path to the HE_scripts housed individually within each sample output directory.

A few places throughout the runHyper script we found it to be helpful to hardcode the path, either to components of the hyper-editing pipeline (lines 48-49) or to the unmapped bam files (lines 79, 82, and 85). These places are also notated with comments throughout the runHyper script itself.

## Additional Processing Steps
#### Step One: Preparing A2G BED files
Since the example we provided with this repository is data from paired-end sequencing, the pipeline will generate two separate output files which must be concatenated together. If single-end sequencing data was used, this isn't necessary. We are only considering the A2G hyper-edited clusters (UE.bed) and sites (ES.bed).

```
PATH1="/full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-1.UE.bed_files"

PATH2="full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-2.UE.bed_files"

cat ${PATH1}/A2G.bed ${PATH2}/A2G.bed > H276_GABA.UE.bed
```
```
PATH1="/full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-1.ES.bed_files"

PATH2="full/path/to/experiment_name/H276_GABA/UEdetect.PE_0.05_0.6_30_0.6_0.1_0.8_0.2/H276_GABA-2.ES.bed_files"

cat ${PATH1}/A2G.bed ${PATH2}/A2G.bed > H276_GABA.ES.bed
```

#### Step Two: Filtering editing sites that may represent common variants or reside within problematic regions of the genome

First, the ES BED file must be converted into a VCF file to input into [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/#table_annovarpl), then we make use of the table_annovar.pl program to annotate common variants aggregated across several databases. We use [Bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) to execute filtering. Blacklisted regions of the genome derived from [ENCODE](https://www.nature.com/articles/s41598-019-45839-z). Here, we additionally filtered out those editing sites which reside in homopolymeric stretches of the genome.

###### BED to VCF
```
ES <- read.delim("H276_GABA.ES.bed", header = FALSE, sep = "\t")
ES[,2] <- ES[,2]+1
Ref <- ES$V6
Ref[Ref == "+"] <- "A"
Ref[Ref == "-"] <- "T"

Alt <- ES$V6
Alt[Alt == "+"] <- "G"
Alt[Alt == "-"] <- "C"

annovarES <- cbind(ES, Ref, Alt)
annovarES <- annovarES[, c(1, 2, 3, 7, 8, 4, 5, 6)]
colnames(annovarES) <- c("Chr", "Start", "Stop", "Ref", "Alt",
                           "Info1", "Info2", "Info3")

write.table(annovarES, file = "anno.H276_GABA.ES.bed,
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
}
```

###### Annotation and filtering

```
IN="/full/path/to/experiment_name/H276_GABA/anno.H276_GABA.ES.bed"
OUT="/full/path/to/experiment_name/H276_GABA.myanno"

table_annovar.pl $IN resources/annovar/humandb/ -buildver hg38 -out $OUT -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920 -operation g,f,f,r,r,f --argument ,,,\'--colsWanted 5\',\'--colsWanted 10&11&12\', -nastring "." --otherinfo --thread 10 --maxGeneThread 10

awk 'BEGIN{OFS=FS="\t"}{if ( ($11=="." || $11=="dbsnp153CommonSNV") && ( $12<0.05 || $12=="AF")) print $0}' H276_GABA.myanno.hg38_multianno.txt > H276_GABA.myanno.hg38_multianno.txt.noCommon.txt

bedtools intersect -a H276_GABA.myanno_multianno.txt.noCommon.txt -b hg38-blacklist.v2_sort.bed, homopolymeric_sites.bed -v > H276_GABA.SNPsremoved_noblacklist.txt

perl format_HE_files.pl H276_GABA.SNPsremoved_noblacklist.txt > H276_GABA.filtered.ES.bed
```

#### Step Three: Merging overlapping hyper-editing clusters and counting number of editing sites per cluster

We merge the hyper-edited clusters (from the UE BED file), so any features with overlapping genomic coordiantes are reported as a single cluster. Then we use the ES BED file to count the number of editing sites per cluster, which is used in the Step Four of filtering. The code below utilizes Bedtools.

```
sort -k1,1 -k2,2n H276_GABA.UE.bed > H276_GABA.UE.bed.sorted

bedtools merge -i H276_GABA.UE.bed.sorted > H276_GABA.UE.bed.sorted.merged

bedtools intersect -a H276_GABA.UE.bed.sorted.merged -b H276_GABA_final.txt -wa -c > H276_GABA.UE.bed.counted.temp
```

#### Step Four: Stretching hyper-editing clusters

On a sample level, the average length (in nucleotides) of each hyper-editing cluster is calculated, in addition to the average number of A2G editing sites within those hyper-editing clusters. This number is then subtracted from the "start" coordinate and added to the "end" coordinate, effectively lengthening the cluster by 2x the average distance between editing sites, in order to account for the possibility that low-sequencing coverage at the edges of the defined cluster borders prevents inclusion of editing events resonably nearby.

```
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

#### Step Four: Final hyper-editing clusters merging

The "stretched" hyper-editing clusters are merged one more time, combining features with overlapping genomic coordinates into one.

```
bedtools merge -i H276_GABA.UE.stretched > H276_GABA.UE.bed.final

```

##### The final outputs from this pipeline are a BED file of hyper-editing clusters and a BED file of individual hyper-editing sites annotated with gene, genic region, and repeat type.
