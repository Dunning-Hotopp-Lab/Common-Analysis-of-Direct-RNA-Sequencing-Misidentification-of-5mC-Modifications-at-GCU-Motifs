# Common Analysis of Direct RNA SequencinG CUrrently leads to Misidentification of 5-Methylcytosine Modifications
Kaylee Watson<br />


## Table of Contents

[• Basecalling ONT Reads with Guppy](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#basecalling-ont-reads-with-guppy)<br />
[• Mapping to Transcriptome References](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#mapping-to-transcriptome-references)<br />
[• Sequencing Stats](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#sequencing-stats)<br />
[• Tombo Modification Detection - Alternative Model](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#tombo-modification-detection---alternative-model)<br />
[• Motif Detection in Top 1000 Modified Positions](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#motif-detection-in-top-1000-modified-positions)<br />
[• Modified Fraction Boxplots](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#modified-fraction-boxplots)<br />
[• Modified Fraction Density Plots - Sindbis virus](https://github.com/kayleewatson/Common-Analysis-of-Direct-RNA-SequencinG-CUrrently-leads-to-Misidentification-of-5-Methylcytosine/edit/main/README.md#modified-fraction-density-plots---sindbis-virus)<br />

## Basecalling ONT Reads with Guppy
#### Command for basecalling RNA using Guppy 6.4.2:
```
FAST5_DIR = path_to_fast5
OUTPUT_DIR = output_directory_path

/usr/local/packages/guppy-6.4.2_gpu/bin/guppy_basecaller -x "cuda:0" --input_path "$FAST5_DIR" --save_path "$OUTPUT_DIR" --config rna_r9.4.1_70bps_hac.cfg --min_qscore 7 --records_per_fastq 10000000 --gpu_runners_per_device 8 --num_callers 1
```


## Mapping to Transcriptome References
Map to transcriptome reference with minimap2 (since Tombo uses transcriptome)
```
REF=path_to_transcriptome_ref
FASTQ=path_to_fastq
SAM_FILE=name_of_SAM_output


BAM_FILE=name_of_BAM_output
samtools view -bhF 2308 $SAM_FILE | samtools sort -o $BAM_FILE
#filter out reads that are not primary alignments
```

## Sequencing Stats

### Determine total reads mapped
```
BAM_FILE=path_to_BAM_file
samtools view $BAM_FILE | wc -l 
```

### Mbp mapped, N50, max read length
Remove any soft-clipped regions for more accurate stats (removes ONT adapter regions)
```
JVARKIT=path_to_jvarkit
BAM_FILE=path_to_BAM_input
BAM_OUT=name_of_BAM_output
/usr/local/packages/jdk-8u171/bin/java  -Xmx10g -jar $JVARKIT/dist/biostar84452.jar -o $BAM_OUT $BAM_FILE
```

Get statistics
```
BAM_OUT=(same as previous step)
FASTQ_OUT=name_of_fastq_output
STATS_FILE=name_of_stats_output

samtools fastq $BAM_OUT > $FASTQ_OUT
# create FASTQ from BAM file to run statistics

seqkit stats -T -a -o $STATS_FILE < $FASTQ_OUT
```

To get percent mapped, pull stats from original fastq file and then use sum_len column from both stats files to get percentage
```
FASTQ=path_to_original_fastq
STATS_FILE=name_of_stats_output

seqkit stats -T -a -o $STATS_FILE < $FASTQ
```

### Percent rRNA
Filter gff file for rRNA regions (column 3)
```
GFF=path_to_gff
GFF_RRNA=name_of_rRNA_gff

awk '$3 ~ /rRNA/ { print }' $GFF > $GFF_RRNA
```
Filter BAM file for reads that mapped to rRNA and get number of reads from that
```
PATH=/usr/local/packages/bedops-2.4.36:"$PATH"
BED=name_of_BED_file
GFF_RRNA=path_to_rRNA_gff
BAM=path_to_BAM_file
BAM_OUT=name_of_BAM_output

gff2bed < $GFF_RRNA > $BED
# Convert gff to BED file

samtools view -L $BED $BAM -o $BAM_OUT
samtools view $BAM_OUT | wc -l
```

## Tombo Modification Detection - Alternative Model

### ONT FAST5 API: Multi-to-Single FAST5
Tombo requires single FAST5 format
```
FAST5_DIR=path_to_fast5_files

OUTPUT_DIR=output_path

multi_to_single_fast5 -i "$FAST5_DIR" -s "$OUTPUT_DIR"
```
### Tombo Reannotate
Reannotate single FAST5 files with basecalls from FASTQ
```
FAST5_DIR=path_to_single_fast5_files
FASTQ_FILE=/path_to_fastq_file

tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir "$FAST5_DIR" --fastq-filenames "$FASTQ_FILE"
```
### Tombo Resquiggle
Resquiggle FAST5 files with raw signal
```
REF_TRANSCRIPT_FNA=path_to_reference

FAST5_DIR=path_to_reannotated_fast5

tombo resquiggle --overwrite --rna "$FAST5_DIR" "$REF_TRANSCRIPT_FNA" --num-most-common-errors 5 --ignore-read-locks
```

### Tombo Detect Modifications
Detect modifications using the Alternative Model
```
FAST5_DIR=path_to_resquiggled_fast5
STATS_PREFIX=prefix_for_stats_files

tombo detect_modifications alternative_model --fast5-basedirs "$FAST5_DIR" --statistics-file-basename "$STATS_PREFIX" --per-read-statistics-basename "$STATS_PREFIX" --alternate-bases 5mC
```

## Motif Detection in Top 1000 Modified Positions

### Tombo Output Significant Regions
```
REF=path_to_reference
FAST5_DIR=path_to_resquiggled_fast5
STATS_FILE=path_to_stats_file (output from detect_modifications)

tombo text_output signif_sequence_context --fast5-basedirs "$FAST5_DIR" --statistics-filename HEK293_SRR13261194.5mC_base_detection.5mC.tombo.stats --num-regions 1000 --num-bases 10
```

#### Change all T's (from Tombo output) to U's in fasta:
```
perl -pe 'tr/tT/uU/ unless(/>/)' < file.fasta > file.RNA.fasta
```

## Modified Fraction Boxplots
#### Tombo - Dampened Fraction of Modified Reads:

```
TOMBO_BIN_DIR = /local/projects-t3/JULIE/packages/miniconda3/envs/epitombo/bin
STATS_FILE = path_to_stats_file
FILE_PREFIX = sample_prefix
REF = path_to_reference_transcriptome

"$TOMBO_BIN_DIR"/tombo text_output browser_files --genome-fasta "$REF" --statistics-filename "$STATS_FILE" --browser-file-basename "$FILE_PREFIX" --file-types dampened_fraction

"$TOMBO_BIN_DIR"/tombo text_output browser_files --genome-fasta "$REF" --statistics-filename "$STATS_FILE" --browser-file-basename "$FILE_PREFIX" --motif-descriptions GCT:2:5mC --file-types dampened_fraction
```

#### Filter by requiggled read depth

Run Tombo to output 'valid coverage' at every position
```
REF=path_to_reference
STATS_FILE=path_to_stats_file
FILE_PREFIX=name_of_output_file

tombo text_output browser_files --genome-fasta "$REF" --statistics-filename "$STATS_FILE" --browser-file-basename "$FILE_PREFIX" --file-types valid_coverage
```
Genes/regions where reads did not map will have blank space after the WIG header. This (and the header for that region) needs to be removed, then create BED file (with all columns) from WIG file:
```
PATH=/usr/local/packages/bedops-2.4.36:"$PATH"

for i in *dampened_fraction_modified_reads.plus.wig; do sed -i '$!N;/=.*\n$/d;P;D' $i; done
# the wig header has "=" in it, so this command finds lines with "=", and if the next line is blank, removes both
for i in *dampened_fraction*.plus.wig; do wig2bed-typical < $i > $i.bed; done

for i in *valid_coverage.plus.wig; do sed -i '$!N;/=.*\n$/d;P;D' $i; done
for i in *valid_coverage.plus.wig; do wig2bed-typical < $i > $i.bed; done
```
VIRUS (low depth samples): Remove lines in valid_coverage that have < 10 in the 5th column (depth column)
```
for i in *.valid_coverage.plus.wig.bed; do awk '{if($5>9) print}' $i > $i.filtered; done
```
ALL OTHERS (B. malayi, D. ananassae, C. albicans, E. coli): Remove lines in valid_coverage that have < 100 in the 5th column
```
for i in *.valid_coverage.plus.wig.bed; do awk '{if($5>99) print}' $i > $i.filtered; done
```
Filter each bed file by valid_coverage file
```
FILT_COV=filtered_valid_coverage_file
DAMP_FRAC=dampened_fraction_file
FILT_FRAC=filtered_output_name

awk 'NR==FNR{a[$2]=1;next}a[$2]' $FILT_COV $DAMP_FRAC > $FILT_FRAC
```

Create a file with Non-GCU motifs
```
for i in *filtered.bed; do awk '{print $1"\t"$2"\t"$3"\t"$5}' $i > $i.4col && mv $i.4col $i; done
# remove the 4th column since it won't match between the 2 files (just corresponds to row number)

ALL_FILE=all_motifs_file
GCU_FILE=GCU_motif_file
PREFIX=output_file_prefix

comm -2 -3 <(sort $ALL_FILE) <(sort $GCU_FILE) > $PREFIX.nonGCU.dampened_fraction_filtered.bed
```

Create BED files formatted for boxplots from filtered files by adding a column for the motif
```
for i in *all.dampened_fraction_filtered.bed; do awk -v motif="All" '{print motif"\t"$4}' $i > $i.formatted; done
for i in *nonGCU.dampened_fraction_filtered.bed; do awk -v motif="Non-GCU" '{print motif"\t"$4}' $i > $i.formatted; done
for i in *_GCU.dampened_fraction_filtered.bed; do awk -v motif="GCU" '{print motif"\t"$4}' $i > $i.formatted; done
```
Add a column with the sample name, then combine all 'final' files in a single file for R:
```
for i in 20190701_Bmalayi*.formatted; do awk -v sample="Bmalayi" '{if($2>0) print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20191024_Calbicans*.formatted; do awk -v sample="Calbicans" '{if($2>0) print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in Dananassae*.formatted; do awk -v sample="Dananassae" '{if($2>0) print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20220727_Ecoli*.formatted; do awk -v sample="Ecoli" '{if($2>0) print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20220512_SINV_IVT*.formatted; do awk -v sample="SINV_IVT" '{if($2>0) print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20181026_JW18*.formatted; do awk -v sample="JW18_SINV" '{if($2>0) print sample"\t"$1"\t"$2}' $i > final.$i; done

cat final* > modified_fractions_all.tsv
```
Check number of positions analyzed for each motif:
```
wc -l final*
```

#### Plot Modified Fractions
Command to run R script on final concatenated file:
```
~/scripts/boxplot_mods_multisample.r modified_fractions_all.tsv
```
Package requirements:
* ggplot2

R script for boxplots:
```
#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please supply a filename", call.=FALSE)
}

cat("Gathering the data...\n")

df <- read.delim(args[1], header=FALSE, sep="\t")
df$V1 <- as.factor(df$V1)

df_nozero <- df[(df$V3 > 0),]
# removes all '0' modified fractions, so only nonzero fractions are plotted

df_nozero$V2 <- factor(df_nozero$V2, levels = c("All", "Non-GCU", "GCU"))

cat("Plotting...\n")

jpeg("plot.jpeg", width=2000, height=550, units="px")
# create a jpeg file with specified size and units

ggplot(df_nozero, aes(x=V2, y=V3, fill=V2)) +
geom_boxplot(outlier.size = 0.3) +
scale_fill_manual(values=c("white","grey70","grey30")) +
facet_grid(~factor(V1, levels=c("Bmalayi", "Dananassae",
    "Calbicans", "Ecoli", "JW18_SINV", "SINV_IVT"))) +
# create separate boxplot for each organism/sample in a specific order
ylab("Fraction of\nReads Modified") +
theme_bw() +
labs(fill="                                                            ") +
# use a large space for the legend title to shift the legend to the bottom right
    theme(
axis.title.y=element_text(size=45),
axis.title.x=element_blank(),
axis.text.x=element_blank(),
axis.ticks.x = element_blank(),
axis.text.y=element_text(size=35),
strip.text.x = element_text(size=25),
legend.text=element_text(size=50),
legend.title=element_text(size=50),
legend.position = "bottom",
legend.direction = "horizontal",
panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
scale_x_discrete(limits = c("All", "Non-GCU", "GCU"))

invisible(dev.off())
```
#### Z-test and Cohen's d


## Modified Fraction Density Plots - Sindbis virus

### WIG to BED

```
WIG_FILE = path_to_wig_file
# this is the output from Tombo
BED_FILE = bed_filename

PATH=/usr/local/packages/bedops-2.4.36:"$PATH"
wig2bed-typical < $WIG_FILE > $BED_FILE 
```

### Remove positions with a fraction of 0.000000 in both IVT and JW18 BED files

```
BED_FILE = path_to_bed_file
NOZERO_BED = output_filename
awk '{if($5 > 0) print}' $BED_FILE > $NOZERO_BED
```

### Remove positions not in both IVT and JW18 BED files

```
NOZERO_BED = path_to_nozero_bed_file
POSITIONS_FILE = output_filename
FILTERED_FILE = final_filtered_filename

awk '{print $1"\t"$2}' $NOZERO_BED > $POSITIONS_FILE
grep -Fwf file1 file2 > common_positions.txt
#file1 and file2 are the POSITIONS_FILE outputs for both JW18 and IVT

awk 'NR==FNR{a[$2];next} $2 in a{print}' common_positions.txt $NOZERO_BED > $FILTERED_FILE
```

### Create file with 3-mer and modified fractions

#### Bash script to create 3-mer file (create_3mer_fraction_file.sh):
```
#!/usr/bin/env bash

#this script requires 3 arguments:
## first argument is the path to the bed file of modified fractions (FILTERED_FILE)
## second argument is the path to the reference fasta
## third argument is the sample name

tmpfile=$(mktemp --suffix=.txt)
tmpfile2=$(mktemp --suffix=.txt)
flanking_tsv=$(mktemp --suffix=.tsv)
tmptsv=$(mktemp --suffix=.tsv)

awk '{print $1">"$2":"$3"<"$4"+"$5}' $1 > $tmpfile

python ~/scripts/flanking_nt_intervals.py $tmpfile > $tmpfile2

bedtools getfasta -fi $2 -bed $tmpfile2 -tab | awk -F ":" 'OFS="\t" {print $1,$2,$3}' | awk -F "-" 'OFS="\t" {print $1,$2,$3,$4}' > $tmptsv

join -j 2 -o 1.5,2.4 <(sort -k 1b,2 $tmpfile2) <(sort -k 1b,2 $tmptsv) | perl -pe 'tr/tT/uU/' | awk -v myvar=$3 '{print $1"\t"$2"\t"myvar}'

rm $tmpfile
rm $tmpfile2
rm $flanking_tsv
rm $tmptsv
```
### Run R script to Plot and Calculate 3-mer Sites

Package requirements:
* ggplot2
* scales

```
#!/usr/bin/env Rscript

library(ggplot2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please supply a filename", call.=FALSE)
}

cat("Reading in the data...")
data1 = read.delim(args[1], header=FALSE, sep="\t")
data2 = read.delim(args[2], header=FALSE, sep="\t")

cat("\nFormatting the data...\n\n")
colnames(data1) <- c("Value", "Identifier", "Sample")
colnames(data2) <- c("Value", "Identifier", "Sample")
data1$Identifier <- as.factor(data1$Identifier)
data2$Identifier <- as.factor(data2$Identifier)
data1$Sample <- as.factor(data1$Sample)
data2$Sample <- as.factor(data2$Sample)

data1_motif_count <- table(data1$Identifier)

cat(as.character(data1$Sample[1]))
print(data1_motif_count)

data2_motif_count <- table(data2$Identifier)

cat(as.character(data2$Sample[1]))
print(data2_motif_count)


cat("\nPlotting...")

jpeg("density_plot.jpeg", width=2000, height=1550, units="px")

ggplot(aes(Value, group = Identifier), data = data1) +
    geom_density(aes(Value, group = Identifier,color = Sample), size=2.5, data = data1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    colour = "black", fill = NA) +
    geom_density(aes(Value, group = Identifier, color = Sample), size=2.5, data = data2) +
    facet_wrap(~Identifier) +
    scale_color_manual(name = "Viral\nRNA", values = c("black", "grey55")) +
    xlab("Methylated Fraction") + theme_linedraw() +
    theme(axis.title=element_text(size=40),
    legend.text=element_text(size=45),
    legend.title=element_text(size=45),
    strip.text = element_text(size=40, face="bold"),
    axis.text=element_text(size=30),
    legend.position="right",
    strip.background = element_rect(color=NA),
    panel.spacing = unit(2, "lines")) +
    coord_cartesian(ylim = c(0, 4)) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,0.25,0.5,0.75,1))

sample_name = as.character(data1$Sample[1])
filename = paste(sample_name, "_RNA_histogram_density.jpeg", sep="")
jpeg(filename, width=2000, height=1550, units="px")

ggplot(aes(Value, group = Identifier), data = data1) +
    geom_histogram(binwidth=0.01, colour = 1, fill = "white") +
    geom_density(aes(Value, group = Identifier,color = Sample), size=2.5, data = data1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    colour = "black", fill = NA) +
    facet_wrap(~Identifier) +
    scale_color_manual(name = "Viral\nRNA", values = c("black", "grey55")) +
    xlab("Methylated Fraction") + theme_linedraw() +
    theme(axis.title=element_text(size=40),
    legend.text=element_text(size=45),
    legend.title=element_text(size=45),
    strip.text = element_text(size=40, face="bold"),
    axis.text=element_text(size=30),
    legend.position="none",
    strip.background = element_rect(color=NA),
    panel.spacing = unit(2, "lines")) +
    coord_cartesian(ylim = c(0, 12)) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,0.25,0.5,0.75,1))

sample_name = as.character(data2$Sample[1])
filename = paste(sample_name, "_RNA_histogram_density.jpeg", sep="")
jpeg(filename, width=2000, height=1550, units="px")

ggplot(aes(Value, group = Identifier), data = data2) +
    geom_histogram(binwidth=0.01, colour = 1, fill = "white") +
    geom_density(aes(Value, group = Identifier,color = Sample), size=2.5, data = data2) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    colour = "black", fill = NA) +
    facet_wrap(~Identifier) +
    scale_color_manual(name = "Viral\nRNA", values = c("black", "grey55")) +
    xlab("Methylated Fraction") + theme_linedraw() +
    theme(axis.title=element_text(size=40),
    legend.text=element_text(size=45),
    legend.title=element_text(size=45),
    strip.text = element_text(size=40, face="bold"),
    axis.text=element_text(size=30),
    legend.position="none",
    strip.background = element_rect(color=NA),
    panel.spacing = unit(2, "lines")) +
    coord_cartesian(ylim = c(0, 12)) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,0.25,0.5,0.75,1))


invisible(dev.off())
cat("\n")
```
