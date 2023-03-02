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

#### Tombo - Dampened Fraction of Modified Reads

#### WIG file to BED file

Genes/regions where reads did not map will have blank space after the WIG header. This (and the header for that region) needs to be removed.

```
WIG_FILE = path_to_wig_file
BED_FILE = name_of_bed_file
# the wig header has "=" in it, so this command finds lines with "=", and if the next line is blank, removes both
sed -i '$!N;/=.*\n$/d;P;D' $WIG_FILE

# change WIG file to BED file format
wig2bed-typical < "$WIG_FILE" > "$BED_FILE"
```

### Modified Fraction Boxplots
File path variables
```
GCU_DIR=/local/projects-t3/RDEPI/GCU_tombo
BMALAYI_ALL=$GCU_DIR/Bmalayi/Bmalayi.dampened_fraction_modified_reads.plus.bed
BMALAYI_NON_GCU=$GCU_DIR/Bmalayi/Bmalayi_HCV.dampened_fraction_modified_reads.5mC.plus.bed
BMALAYI_GCU=$GCU_DIR/Bmalayi/Bmalayi_GCU.dampened_fraction_modified_reads.5mC.plus.bed

DANA_ALL=$GCU_DIR/Dananassae/Dananassae.dampened_fraction_modified_reads.plus.bed
DANA_NON_GCU=$GCU_DIR/Dananassae/Dananassae_HCV.dampened_fraction_modified_reads.5mC.plus.bed
DANA_GCU=$GCU_DIR/Dananassae/Dananassae_GCU.dampened_fraction_modified_reads.5mC.plus.bed

CALBICANS_ALL=$GCU_DIR/Calbicans/Calbicans.dampened_fraction_modified_reads.plus.bed
CALBICANS_NON_GCU=$GCU_DIR/Calbicans/Calbicans_HCV.dampened_fraction_modified_reads.5mC.plus.bed
CALBICANS_GCU=$GCU_DIR/Calbicans/Calbicans_GCU.dampened_fraction_modified_reads.5mC.plus.bed

ECOLI_ALL=$GCU_DIR/K12_LB/K12_LB.dampened_fraction_modified_reads.plus.bed
ECOLI_NON_GCU=$GCU_DIR/K12_LB/K12_LB_HCV.dampened_fraction_modified_reads.5mC.plus.bed
ECOLI_GCU=$GCU_DIR/K12_LB/K12_LB_GCU.dampened_fraction_modified_reads.5mC.plus.bed

SINV_ALL=$GCU_DIR/SINV/JW18_SINV.dampened_fraction_modified_reads.plus.bed
SINV_NON_GCU=$GCU_DIR/SINV/JW18_HCV.dampened_fraction_modified_reads.5mC.plus.bed
SINV_GCU=$GCU_DIR/SINV/JW18_GCU.dampened_fraction_modified_reads.5mC.plus.bed

IVT_SINV_ALL=$GCU_DIR/SINV/SINV_IVT.dampened_fraction_modified_reads.plus.bed
IVT_SINV_NON_GCU=$GCU_DIR/SINV/IVT_HCV.dampened_fraction_modified_reads.5mC.plus.bed
IVT_SINV_GCU=$GCU_DIR/SINV/IVT_GCU.dampened_fraction_modified_reads.5mC.plus.bed
```

Command to run R script:
```
~/scripts/boxplot_mods_wrap.r $BMALAYI_ALL $BMALAYI_NON_GCU $BMALAYI_GCU $DANA_ALL $DANA_NON_GCU $DANA_GCU $CALBICANS_ALL $CALBICANS_NON_GCU $CALBICANS_GCU $ECOLI_ALL $ECOLI_NON_GCU $ECOLI_GCU $SINV_ALL $SINV_NON_GCU $SINV_GCU $IVT_SINV_ALL $IVT_SINV_NON_GCU $IVT_SINV_GCU
```


### Modified Fraction Density Plots - Sindbis virus
