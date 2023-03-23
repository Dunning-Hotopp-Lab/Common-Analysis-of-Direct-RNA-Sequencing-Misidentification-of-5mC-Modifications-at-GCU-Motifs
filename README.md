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

## Sequencing Stats

## Tombo Modification Detection - Alternative Model

## Motif Detection in Top 1000 Modified Positions

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

"$TOMBO_BIN_DIR"/tombo text_output browser_files --genome-fasta "$REF" --statistics-filename "$STATS_FILE" --browser-file-basename "$FILE_PREFIX" --motif-descriptions HCV:2:5mC --file-types dampened_fraction
```

#### Format WIG file

Genes/regions where reads did not map will have blank space after the WIG header. This (and the header for that region) needs to be removed:

```
PATH=/usr/local/packages/bedops-2.4.36:"$PATH"
WIG_FILE = wig_file

# the wig header has "=" in it, so this command finds lines with "=", and if the next line is blank, removes both
sed -i '$!N;/=.*\n$/d;P;D' $WIG_FILE
```

Change the WIG files to BED file format, print only the modified fraction, and add a column for the 3-mer:
```
for i in *GCU.dampened_fraction_modified_reads.5mC.plus.wig; do wig2bed-typical < $i | awk -v motif="GCU" '{print motif"\t"$5}' > $i.bed; done
for i in *HCV.dampened_fraction_modified_reads.5mC.plus.wig; do wig2bed-typical < $i | awk -v motif="Non-GCU" '{print motif"\t"$5}' > $i.bed; done
for i in *dampened_fraction_modified_reads.plus.wig; do wig2bed-typical < $i | awk -v motif="All" '{print motif"\t"$5}' > $i.bed; done
```
Add a column with the sample name, then combine all 'final' files in a single file for R:
```
for i in 20190701_Bmalayi*.bed; do awk -v sample="Bmalayi" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20191024_Calbicans*.bed; do awk -v sample="Calbicans" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in Dananassae*.bed; do awk -v sample="Dananassae" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20220727_Ecoli*.bed; do awk -v sample="Ecoli" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20220512_SINV_IVT*.bed; do awk -v sample="SINV_IVT" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20181026_JW18*.bed; do awk -v sample="JW18_SINV" '{print sample"\t"$1"\t"$2}' $i > final.$i; done

cat final* > modified_fractions_all.tsv
```
#### Filter by requiggled read depth of 10

Run Tombo to output 'valid coverage' at every position
```
REF=path_to_reference
STATS_FILE=path_to_stats_file
FILE_PREFIX=name_of_output_file

tombo text_output browser_files --genome-fasta "$REF" --statistics-filename "$STATS_FILE" --browser-file-basename "$FILE_PREFIX" --file-types valid_coverage
```
Remove blank lines and change wig file to bed file format
```
for i in *valid_coverage.plus.wig; do sed -i '$!N;/=.*\n$/d;P;D' $i; done
for i in *valid_coverage.plus.wig; do wig2bed-typical < $i > $i.bed; done
```
Create BED file (with all columns) from dampened_fraction WIG file:
```
for i in *dampened_fraction_modified_reads.plus.wig; do sed -i '$!N;/=.*\n$/d;P;D' $i; done
for i in *dampened_fraction*.plus.wig; do wig2bed-typical < $i > $i.bed; done
```
Remove lines in valid_coverage that have < 10 in the 5th column (depth column)
```
for i in *.valid_coverage.plus.wig.bed; do awk '{if($5>9) print}' > $i.filtered; done

```
Filter each bed file by valid_coverage file
```
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
#### Wilcoxon–Mann–Whitney test


## Modified Fraction Density Plots - Sindbis virus
