# Transposon-mobility
Collection of functions, scripts, and notebooks used in genome-wide analysis of heat stress-stimulated transposon mobility.


## Ordering and description of scripts

[transposonmobility](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/transposonmobility.py) : A library of python functions used in subsequent analysis.


[Calculate_sample_depth](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/Calculate_sample_depth.ipynb) : Formats *bwa* and *samtools* commands for calculating genome-wide read depth across samples.


[bwa_XL280np](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/bwa_XL280np.sh) : Calls *bwa* and *samtool* commands to align sequenced reads and calculate depth. This script is generated from the *Calculate_sample_depth notebook* and will need to be regenerated after gather fastq files from the *SRA* (see [FASTQ](https://github.com/magwenelab/Transposon-mobility/tree/main/FASTQ))


[Plot_depth](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/Plot_depth.ipynb) : Generates diagnostic plots of genome-wdie sequencing depth per sample.


[Checking_XL280_NP_GFF](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/Checking_XL280_NP_GFF.ipynb) : Checks a few of the gene coding sequences of the XL280 nanopore GFF given the example JEC21 gene sequences.


[Distribution_of_tandom_genes_XL280_NP](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/Distribution_of_tandom_genes_XL280_NP.ipynb) : Calculates the distributions of gene orientations within the XL280 nanopore GFF.


[LocaTE_BLAT](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/LocaTE_BLAT.ipynb) : Formats and submits *blat* commands per sample and locates candidate transposable element insertion sites.


[Plot_TE](https://github.com/magwenelab/Transposon-mobility/blob/main/SCRIPTS/Plot_TE.ipynb) : Plots results of *blat* submissions and searches for transposable element insertion and deletion sites. 


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

### Python packages
    
    os 
    glob
    gzip
    numpy (v 1.19.2)
    pandas (v 1.1.3)
    shutil (v 1.0.0)
    biopython (v 1.78)
    matplotlib (v 3.3.1)
    subprocess