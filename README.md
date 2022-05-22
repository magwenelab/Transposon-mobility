# Transposon-mobility
Collection of functions, scripts, and notebooks used in genome-wide analysis of heat stress-stimulated transposon mobility.

## Dependicies

### BLAT
Install the [Blast-Like Alignment Tool](https://genome.cshlp.org/content/12/4/656.long) via conda:

    conda install -c bioconda blat
    
### Samtools
Install [samtools](www.htslib.org) via conda:

    conda install -c bioconda samtools=1.3.1 

### BWA
Install [BWA](http://bio-bwa.sourceforge.net/) via conda:

    conda install -c bioconda bwa=1.3.1

### Needed python packages
    
    os 
    glob
    gzip
    numpy (v 1.19.2)
    pandas (v 1.1.3)
    shutil (v 1.0.0)
    biopython (v 1.78)
    matplotlib (v 3.3.1)
    subprocess