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

### Modified Fraction Density Plots - Sindbis virus
