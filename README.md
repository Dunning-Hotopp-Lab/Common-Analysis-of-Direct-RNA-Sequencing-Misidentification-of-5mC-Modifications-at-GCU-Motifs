# GCU Paper
Kaylee Watson<br />
01 March 2023

## Table of Contents

[• Basecalling ONT Reads with Guppy](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#basecalling-ont-reads-with-guppy)<br />
[• Mapping to Transcriptome References](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#mapping-to-transcriptome-references)<br />
[• Sequencing Stats](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#sequencing-stats)<br />
[• Tombo Modification Detection - Alternative Model](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#tombo-modification-detection---alternative-model)<br />
[• Motif Detection in Top 1000 Modified Positions](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#motif-detection-in-top-1000-modified-positions)<br />
[• Modified Fraction Boxplots](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#modified-fraction-boxplots)<br />
[• Modified Fraction Density Plots - Sindbis virus](https://github.com/kayleewatson/GCU-Paper/new/main?readme=1#modified-fraction-density-plots---sindbis-virus)<br />

### Basecalling ONT Reads with Guppy

### Mapping to Transcriptome References

### Sequencing Stats

### Tombo Modification Detection - Alternative Model

### Motif Detection in Top 1000 Modified Positions





### Modified Fraction Boxplots
#### Tombo - Dampened Fraction of Modified Reads

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

Genes/regions where reads did not map will have blank space after the WIG header. This (and the header for that region) needs to be removed.

```
WIG_FILE = wig_file

# the wig header has "=" in it, so this command finds lines with "=", and if the next line is blank, removes both
sed -i '$!N;/=.*\n$/d;P;D' $WIG_FILE
```

Change the WIG files to BED file format, print only the modified fraction, and add a column for the 3-mer
```
for i in *GCU.dampened_fraction_modified_reads.5mC.plus.wig; do wig2bed-typical < $i | awk -v motif="GCU" '{print motif"\t"$5}' > $i.bed; done
for i in *HCV.dampened_fraction_modified_reads.5mC.plus.wig; do wig2bed-typical < $i | awk -v motif="Non-GCU" '{print motif"\t"$5}' > $i.bed; done
for i in *dampened_fraction_modified_reads.plus.wig; do wig2bed-typical < $i | awk -v motif="All" '{print motif"\t"$5}' > $i.bed; done
```
Add a column with the sample name, then combine all 'final' files in a sinlge file for R
```
for i in 20190701_Bmalayi*.bed; do awk -v sample="Bmalayi" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in Calbicans*.bed; do awk -v sample="Calbicans" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in Dananassae*.bed; do awk -v sample="Dananassae" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20220727_Ecoli*.bed; do awk -v sample="Ecoli" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20220512_SINV_IVT*.bed; do awk -v sample="SINV_IVT" '{print sample"\t"$1"\t"$2}' $i > final.$i; done
for i in 20181026_JW18*.bed; do awk -v sample="JW18_SINV" '{print sample"\t"$1"\t"$2}' $i > final.$i; done

cat final* > modified_fractions_all.tsv
```

#### Plot Modified Fractions
Command to run R script on final concatenated file:
```
~/scripts/boxplot_mods_multisample.r modified_fractions_all.tsv
```


### Modified Fraction Density Plots - Sindbis virus
