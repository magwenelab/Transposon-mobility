## Functions for TE movement analysis
## Author: Cullen Roth, PhD
## Email: cullenroth@gmail.com

## Note for user:
## In addition to having the needed python packages below ...
## be sure to also install via conda "blat".
## For example:
## 
##     conda install -c bioconda blat
##

## Load in pandas, numpy and other usefull mods
import pandas as pd, numpy as np, gzip, subprocess, os, shutil, glob

## Load in the bio package
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## Load in matplotlib
from matplotlib import pyplot as plt

## Set paths of used variables
## Set the path to the plat execuatlable
blatpath = 'blat' #'/home/croth/./blat'

## Set window size (in bp) used in analysis
mysize = 10000

## Set the path to the centromere
thecentpath = '../DATA/XL280NP_Final_CENs.tsv.gz'

## Write ftn for loading in centromeres
def loadcents(centpath = thecentpath):
    
    ## Bring in centromere dataframe
    centromeres = pd.read_csv(centpath,sep='\t')

    ## add contig
    centromeres['Contig'] = [c.split('_')[-1] for c in centromeres['Sequence Name'].tolist()]

    ## Return
    return centromeres

## Set the chromosome map path
thechrompath = '../DATA/chromosome.map.tsv.gz'

## Write the ftn for loading in chromosome path
def loadmap(chrompath = thechrompath):
    return pd.read_csv(chrompath,sep='\t')

## set the path to th gff 
thegffpath = '../DATA/XL280_NP_NUCLEAR_FINAL_gff.csv.gz'

## Bring gff into python
def loadgff(gff_path=thegffpath):
    gff = pd.read_csv(gff_path)
    gff['chrom'] = [a.split('_')[-1] for a in gff.chrom.values]

    ## Sort values
    gff.sort_values(['chrom','Zstart'],inplace=True)

    ## Reset index
    gff.reset_index(drop=True,inplace=True)

    return gff

## Set sample information, passaging, and names as defined by AG
sampledata = [(23,0,4),(21,0,1),(1,30,1),(2,30,1),
              (3,30,1),(4,30,1),(5,30,1),(22,0,2),
              (6,30,2),(7,30,2),(8,30,2),(9,30,2),
              (10,30,2),(11,37,1),(12,37,1),(13,37,1),
              (14,37,1),(15,37,1),(16,37,2),(17,37,2),
              (18,37,2),(19,37,2),(20,37,2)]

## Make into a data frame and set columns 
sampleinfo = pd.DataFrame(sampledata,columns=['Sample','Passaged','Isolate'])

## Set sample manuscript name
sampleinfo['Name'] = ['XL280.3','XL280.1','30-01','Duplicate-30-01','30-02','30-16','Duplicate-30-16','XL280.2','30-17','30-03','30-18','30-19','30-20','37-01','37-02','37-16','37-17','37-03','37-18','37-19','37-20','37-21','37-22']

## Set the color map to be used
colormap = ['#006BA4', '#FF800E','#ABABAB', '#595959','#5F9ED1', 
            '#C85200','#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']

## Define columns used in psl dataframe
pslcolumns = ['match','mismatch','rep_match', 'Ns','Q_gap_count','Q_gap_bases',
              'T_gap_count','T_gap_bases','T_strand','Q_name', 'Q_size', 'Q_start', 
              'Q_end', 'T_name','T_size', 'T_start', 'T_end',
              'block_count', 'block_sizes', 'qStarts', 'tStarts']

## Write functions for parseing and printing blat files
def parse_psl(path, sr = 5, hd = None, rs=6, psl_cols = pslcolumns):
    """Brings in a PSL file from BLAT."""
    ## Load in psl via pandas
    psl_df = pd.read_csv(path,sep='\t',skiprows=sr,header=hd)

    ## define columns
    psl_df.columns = psl_cols

    ## Calculate Q fraction
    psl_df['Q_frac'] = np.round(psl_df['match']/psl_df['Q_size'],rs)
    
    ## Sort dataframe by Q fraction and reset index
    psl_df.sort_values('Q_frac',ascending=False,inplace=True)
    psl_df.reset_index(drop=True,inplace=True)

    ## Return psl dataframe
    return psl_df

def blatem(database, query, pslout, blatloc = blatpath):
    """Submit BLAT command to command line."""
    ## Write a blat command
    blatcommand = '%s %s %s %s'%(blatloc,query,database,pslout)

    ## make into a process 
    process = subprocess.Popen(blatcommand.split(), stdout=subprocess.PIPE)

    ## Communicate that process
    output, error = process.communicate()

    ## Return output and error
    return output,error

def read_fastq(fname):
    """Parses a fastq file."""
    return list(zip(*[iter(''.join([f for f in gzip.open(fname,'rt').read()]).split('\n'))]*4))

def fastout(fileout, records, ftype = 'fasta'):
    """Wrtie out fasta records."""
    ## write to output
    with open(fileout, "w") as output_handle:
        
        ## Write the recors to otuput handel 
        SeqIO.write(records, output_handle, ftype)

        ## Close the handel 
        output_handle.close()
    pass                        

def fastq_to_fasta(fname, fout, ix = None, error1 = "ERROR: expecting a list of file names"):
    """Take fastq sequences and write them to fasta file."""
    ## If fname is a str load the fastq files
    if (type(fname) != 'list') and (type(fname) == str):
    
        ## Load the file 
        k = read_fastq(fname)                        
    else:

        ## otherwaise take the records we all ready have
        assert type(fname) == list, error1
        k = fname                                    
     
    ## If the index isn't defined, define it.  
    if ix is None:
        ix = np.arange(len(k))                       
    
    ## Make each record a sequence record
    recs = [SeqRecord(Seq(r[1]),id='%s'%ix[i]) for i,r in enumerate(k)]                 
    
    ## Write the fasta records to file
    fastout(fout,recs)
    
    ## return the output file name
    return fout

def makedir(path):
    """Makes a directory on path"""
    ## IF the path exists
    if os.path.exists(path):
    
        ## Remove the directory and every in it
        shutil.rmtree(path, ignore_errors=True)

        ## remake the 
        os.makedirs(path, exist_ok=True) 
    
    else: ## Make the directory otherwise
        os.makedirs(path, exist_ok=True)                                 
    pass

def TEanalysis(psl_path):
    """Parses PSL file given path and returns maximum alignment."""
    ## Parse the psl given from input path
    psl = parse_psl(psl_path) 
    
    ## Gather the unique quary hit names
    hits = psl.Q_name.unique() 
    
    ## Take the maximum hits, and hit names
    return psl.loc[np.array([psl[(psl.Q_name == h)]['Q_frac'].idxmax() for h in hits]),:], hits

def blatTE(fname1, fname2, reference, tepath, psltemp='./psltmp', blatdir = '../BLAT', verbose=False):
    """Given a set of paths for paired-end reads (in fastq.gz format), blat reads to a TE and reference path."""
    ## Remove previous runs and remake temp directory
    makedir(psltemp)
    
    ## Write file names and paths, for read one 
    file1, file2 = fname1.split('/')[-1], fname2.split('/')[-1]

    ## make fasta output files
    fout1, fout2 = '%s/%s.fasta'%(psltemp,file1), '%s/%s.fasta'%(psltemp,file2)

    ## Write psl output names
    psl_out1, psl_out2 = '%s/%s.psl'%(psltemp,file1), '%s/%s.psl'%(psltemp,file2)
    
    ## Make compliment read file names
    R1_compliment_out = fout2.split('.fastq.gz')[0]+'_R1_compliment.fasta'
    R2_compliment_out = fout1.split('.fastq.gz')[0]+'_R2_compliment.fasta'
    
    ## Make compliment psl path out
    R1_compliment_out_psl = R1_compliment_out.split('.fasta')[0]+'.psl'
    R2_compliment_out_psl = R2_compliment_out.split('.fasta')[0]+'.psl'

    ## split out the TE name via the path
    theTE = tepath.split('.fasta')[0].split('/')[-1]

    ## Gather sample and lane info
    sample, lane = file1.split('_')[1], file1.split('_')[2][-1]

    ## if in debug mode print to screen
    if verbose:
        print("Locating %s in sample %s from lane %s."%(theTE,sample,lane))
    
    ## Convert fastq to fasta sequences
    fastq_to_fasta(fname1,fout1) 
    fastq_to_fasta(fname2,fout2)
    
    ## Blat reads to the TE 
    blatem(fout1, tepath, psl_out1) 
    blatem(fout2, tepath, psl_out2)
    
    ## Analysis of read one to the TE 
    R1_psl, R1_TE_hits = TEanalysis(psl_out1)
    
    ## Analysis fo read two to the TE 
    R2_psl, R2_TE_hits = TEanalysis(psl_out2)
    
    ## Gather unique read hits to the TE 
    RS = np.unique(np.concatenate([R1_TE_hits, R2_TE_hits]))
    
    ## Take those supporting reads: i.e. both R1 and R2 map to the TE 
    support = [i for i in RS if (i in R1_TE_hits) and (i in R2_TE_hits)]

    ## Complement analysis gather the index of read  
    R1_compliment_ix = [i for i in R1_TE_hits if i not in support]
    R1_compliment = [r for i,r in enumerate(read_fastq(fname2)) if i in R1_compliment_ix]

    ## Check our work
    assert len(R1_compliment) == len(R1_compliment_ix), "ERROR: We are missing reads within compliment one!"

    R2_compliment_ix = [i for i in R2_TE_hits if i not in support]
    R2_compliment = [r for i,r in enumerate(read_fastq(fname1)) if i in R2_compliment_ix]

    ## Check our wrok
    assert len(R2_compliment) == len(R2_compliment_ix), "ERROR: We are missing reads within compliment two!"
    
    ## Check our work for compliment two
    for r in R2_compliment_ix:
        assert r not in R1_compliment_ix, "ERROR: A read in R2 appears in both the anchor and supportive sets!"

    ## Check our work for compliment one
    for r in R1_compliment_ix:
        assert r not in R2_compliment_ix, "ERROR: A read in R1 appears in both the anchor and supportive sets!"

    ## Convert complement fastq to fasta
    fastq_to_fasta(R1_compliment, R1_compliment_out, R1_compliment_ix)
    fastq_to_fasta(R2_compliment, R2_compliment_out, R2_compliment_ix)
    
    ## Blat R1 and R2 complement fasta to XL280 reference
    blatem(R1_compliment_out, reference, R1_compliment_out_psl)
    blatem(R2_compliment_out, reference, R2_compliment_out_psl)
    
    ## Bring in compliment alignments
    R2_compliment_psl = parse_psl(R2_compliment_out_psl)
    R1_compliment_psl = parse_psl(R1_compliment_out_psl)
    
    ## Check work
    for a in R2_compliment_psl.Q_name.unique():
        assert a in R2_compliment_ix
        
    ## Gather paired info from TE hit for R1  
    R1_TE = R1_psl[(R1_psl.Q_name.isin(R1_compliment_psl.Q_name))][['Q_name','T_start','T_end','T_name','Q_frac']]
    R1_TE.columns = ['Q_name'] + ['TE_%s'%r for r in R1_TE.columns[1:]]
    R1_res = R1_compliment_psl.merge(R1_TE)
    R1_res['Read'] = 1

    ## Gather paired info from TE hit for R2
    R2_TE = R2_psl[(R2_psl.Q_name.isin(R2_compliment_psl.Q_name))][['Q_name','T_start','T_end','T_name','Q_frac']]
    R2_TE.columns = ['Q_name'] + ['TE_%s'%r for r in R2_TE.columns[1:]]
    R2_res = R2_compliment_psl.merge(R2_TE)
    R2_res['Read'] = 2
    
    ## Check our work
    for a in R1_res.Q_name.unique():
        if a in R2_res.Q_name.unique():
            print('Error in file %s'%fname1)
    
    ## Concatonate hits
    TE_res = pd.concat([R1_res,R2_res])
    TE_res.reset_index(drop=True,inplace=True)
    
    ## Calculate number of hits per read
    Nhits_R = TE_res.groupby('Q_name').count()[['match']].reset_index().copy()
    Nhits_R.columns = ['Q_name','Nhits']
    
    ## Save final results
    ## Get unique alignment hits acorss genome and paired 
    max_TE_res = TE_res[(TE_res.Q_name.isin(Nhits_R[(Nhits_R.Nhits==1)].Q_name))].copy()
    
    ## Add info and save out csv file. 
    max_TE_res['Sample'] = sample
    max_TE_res['Lane'] = lane
    
    ## Make the ouput blat dir
    os.makedirs('%s/%s'%(blatdir, theTE), exist_ok=True)
    
    ## Save out the blat results
    max_TE_res.to_csv('%s/%s/BLAT_%s_%s_%s.csv.gz'%(blatdir, theTE, theTE, sample, lane))
    
    ## Clean up
    shutil.rmtree(psltemp, ignore_errors=True)
    pass

## Ftns for visulization of resluts
## Write a ftn for returing a TE hit in a psl dataframe
def return_hit(hits,s,x,y,c):
    """Returns a given TE hit (in hits) on a chromosome (c) between x and y (bp) for sample s."""
    return hits[(hits.Sample==s) & (hits.T_start>=x) & (hits.T_end<=y) & (hits.Chrom==c)]

## Ftn for gathering genes (from a gff) near a point on a chromosome
def gather_genes(gff,chm,point):
    """Gathers the genes from a gff within a given set of chromosome coordiantes (bp)."""
    
    ## Initate list
    cg = []
    
    ## For the positive or negative strands 
    for s in ['-','+']:
        
        ## Take all the genes on this strand  with a start less than the point of interest
        temp = gff[(gff.chrom==chrom) & (gff.strand==s) & (gff.Zstart<=point) & (gff.type=='gene')].sort_values('Zstart').ID.values
        
        ## If there are genes in this list append the last (the closest) to list
        if len(temp)>0:
            cg.append(temp[-1])
        
        ## Take all the genes on this strand, with an end less than this point
        temp = gff[(gff.chrom==chrom) & (gff.strand==s) & (gff.Zend<=point) & (gff.type=='gene')].sort_values('Zstart').ID.values
        
        ## If there are genes in this list append the last (the closest) to list
        if len(temp)>0:
            cg.append(temp[-1])
            
        ## Take all the genes on this strand, with a start position greater than this point
        temp = gff[(gff.chrom==chrom) & (gff.strand==s) & (gff.Zstart>=point) & (gff.type=='gene')].sort_values('Zstart').ID.values
        
        ## If there are any such genes append the sorted first to list
        if len(temp)>0:
            cg.append(temp[0])
            
        ## Take all the genes on this strand, with an end position greater than this point
        temp = gff[(gff.chrom==chrom) & (gff.strand==s) & (gff.Zend>=point) & (gff.type=='gene')].sort_values('Zstart').ID.values
        
        ## If there are any such genes append the sorted first to list
        if len(temp)>0:
            cg.append(temp[0])
      
    ## Gather the unique gene names
    cg = np.unique(cg)
      
    ## Slice them from the gff for the negatively oriented set
    gffn = gff[(gff.ID.isin(cg)) & (gff.chrom==chrom) & (gff.strand=='-')]

    ## Slice them from the gff for the positively oriented set
    gffp = gff[(gff.ID.isin(cg)) & (gff.chrom==chrom) & (gff.strand=='+')]
    
    ## return them
    return gffn, gffp

## A ftn for ploting gene models
def plot_genes(tempgff,y1,y2,color,ymod=0,fs=10):
    
    ## Gather genes
    uniquegenes = tempgff[(tempgff.type=='gene')].ID.unique()
    
    ## How many?
    ngenes = len(uniquegenes)
    
    ## Check our work
    assert (len(ngenes) > 0), "ERROR: No genes where given to plot_genes function!"
    
    ## Set the unique y values
    y = np.linspace(y1,y2,ngenes)
    
    ## plot the gene mods
    for i,g in enumerate(uniquegenes):
        ## Gather the gene body
        genebody = tempgff[(tempgff.ID==g) & (tempgff.Type=='gene')][['Zstart','Zend']].values
        
        ## plot the horizontal gene line
        plt.hlines(y[i],np.min(genebody),np.max(genebody),color=color,linewidth=0.5)
    
        ## Gather the CDS of the gene
        cds = tempgff[(tempgff.ID==g) & (tempgff.type=='CDS')][['Zstart','Zend']].values
        
        ## Iterate thru the cds and plot the gene bodies
        [plt.hlines(y[i],l,r,color=color,linewidth=5) for (l,r) in cds]
        
        ## annotate the gene name
        plt.text(np.mean(genebody),y[i]+ymod,s=g,color=color,fontsize=fs, va='center', ha='center')
    pass
    
## A ftn for ploting TE hits
def plot_hit(hits,si,x1,x2,chrom,clen,gff,ax,fs=10,posy=(1700,2300),negy=(300,1500)):
    """Plots the locations of anchor read hits and nerby genes."""
    ## Gather the hit
    sh = return_hit(hits,si,x1,x2,c)
    
    ## Set to the passed axis
    plt.sca(ax)
    
    ## Plot the hit
    plt.plot(sh[['T_start','T_end']].T.mean(),sh[['TE_T_start','TE_T_end']].T.mean(),'kx',alpha=0.25);
    
    ## Plot the mean of the hit
    test = sh[['T_start','T_end']].T.mean().mean()
    
    ## Annotate the number of reads that support the hit
    plt.text(x=0.1,y=0.9, s=str(sh.shape[0]), ha='center', va='center', transform=ax.transAxes)
    
    ## Gather the genes nearest the point in test
    tcn,tcp = gather_genes(gff,chrom,test)
    
    ## Plot the positive genes
    plot_genes(tcp,*posy,'blue', ymod=400, fs = fs)
    
    ## plot the negative genes
    plot_genes(tcn,*negy,'red', ymod=-300, fs = fs)
    pass

## Write a ftn for loading in TE hit dataframe
def loadhits(wildpath,chrmap=None,verbose=True,chrnames=['T_name','T_size','Chrom'],merge_chrommap=True):
    """
    Loads and concatonates TE psl paths across samples.
    This function assumes wildpath is a path with a wild card *
    """
    ## Load and concatonate the dataframes for the paths in wildpath 
    tehits = pd.concat([pd.read_csv(k,index_col=0) for k in sorted(glob.glob(wildpath))])

    ## if in verbose mode, print the number of unique samples to screen
    if verbose:
        print('We have detected data for %s samples'%tehits.Sample.unique().shape[0])


    ## Merge the hits with the chrommap
    if merge_chrommap:
        
        ## Merge the chormosome map
        if chrmap is None:
            chrmap = loadmap()

        tehits = tehits.merge(chrmap[chrnames])

    ## Correct the sample name
    tehits['Sample'] = [int(s[1:]) if s[0] == 'S' else s for s in tehits.Sample.tolist()]

    ## Return the hits
    return tehits

## Write a ftn to make windows
def windows(clen,wsize):
    return np.arange(0,clen+2*wsize,wsize)

## Write a ftn for counting hits via genomic windows 
def tedfall(hits=None,saminfo=None,wildpath=None,chrmap=None,wsize=mysize,verbose=False,offset=mysize/2):
    """
    Given a windowsize counts the number of anchoring reads in that loci.
    """
    ## Gather TE hits
    if hits is None:
        assert (wildpath is not None), "If a dataframe with TE locations (hits) is not given the wildpath must be defined"

        ## load in hits
        hits = loadhits(wildpath, verbose=False)

    ## Gather sample info
    if saminfo is None:
        saminfo = sampleinfo

    ## bring in the windows
    if chrmap is None:
        chrmap = loadmap()

    ## Gather contigs
    contigs = chrmap.Chrom

    ## Iterate thru samples
    ## Initilize dataframe
    dfs = []

    ## Iterate thru the sample info dataframe
    for k,j in saminfo.iterrows():

        ## print the sample
        if verbose:
            print(j.Sample,type(j.Sample))

        ## Take hits for just this sample
        sh = hits[(hits.Sample==j.Sample)]

        assert (sh.shape[0] > 0), "ERROR: Data for this sample is empty!"

        ## Iterate thru congits
        for i,cont in enumerate(contigs):
            
            ## take hits for this chromosome, sorted by left most positons
            shc = sh[(sh.Chrom==cont)].sort_values('T_start')
            
            ## Gather the cumulative positons and chromosome length
            cumpos = chrmap[(chrmap.Chrom==cont)]['Cumpos'].min()
            clen = chrmap[(chrmap.Chrom==cont)]['T_size'].min()
    
            ## Take the te coordinates of anchor reads
            x = shc[['T_start','T_end']].T.mean(axis=0).values

            ## Make windows
            wss = windows(clen,wsize)
    
            ## count the anchor reads
            y = np.array([np.sum((x<=wss[i+1]) & (x>=wss[i])) for i in range(len(wss)-1)])
            
            ## Make a dataframe of hits
            df = pd.DataFrame((wss[:-1],y,np.repeat(cont,len(wss)-1))).T
            
            ## Set column names
            df.columns = ['Pos','Reads','Chrom']

            ## Add cumlative position
            #df['Cumpos'] = df.Pos.values + cumpos
            
            ## Aadd other sample column info
            for i,l in j.iteritems():
                df[i] = l
            
            ## Add to data frame
            dfs.append(df)

    ## concatonate dataframe
    tedfall = pd.concat(dfs)

    ## Adjust positional offset
    tedfall['Pos'] = tedfall['Pos'].values + offset

    ## Check our work
    assert (tedfall.Sample.unique().shape[0] == sampleinfo.shape[0]), "ERROR: We are missing sample data!"

    ## Return dataframe
    return tedfall

## Slice out by chromosome
def bychrom(df,c,n='Chrom'):
    return df[(df[n]==c)]

## Write ftn for plotting TE hits, genome-wide
def manhattan_plot(hit,x='Pos',y='Reads',cumx='Cumpos',colors = colormap,fs=10,figsize=(10,3),
                   savepath=None,xlims=None,ylims=None,chrmap = None,ax=None,
                   xlabel='Chromosome',ylabel='Number of Reads',linewidth=0.5,rast=True):
    """From the dataframe hit, plots the nubmer of TE reads (labeled in y) along the chromosome positions (stored in x)."""
    
    ## call an axis if it is given
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=figsize)
        fig.set_facecolor('w')

    ## Load chrommap
    if chrmap is None:
        chrmap = loadmap()
       
    ## Iterate thru the chromosomes
    for i,r in chrmap.iterrows():
    
        ## Slice by chromosome
        temp = bychrom(hit,r.Chrom)

        ## Calc the cumlative pos
        cumpos = r[cumx]

        ## Cather chrom color
        chromcolor = colors[int(i%len(colors))]

        ## Plot the chromosome body
        plt.hlines(0, np.min(temp[x])+cumpos, np.max(temp[x])+cumpos, 
                   linewidth=linewidth, rasterized=rast, color=chromcolor)

        ## Plot the hits
        plt.vlines(temp[x].values+cumpos, np.zeros(len(temp)), temp[y].values,
                   linewidth=linewidth, rasterized=rast, color=chromcolor)

    ## Adjust xticks and set lables
    plt.xticks(chrmap.Midpts.values,chrmap.index.values+1,fontsize=fs)
    plt.xlabel(xlabel,fontsize=fs)
    
    ## Adjust yticks
    plt.yticks(fontsize=fs)
    plt.ylabel(ylabel,fontsize=fs)
    
    ## Set ylims
    if (ylims is not None):
        plt.ylim(ylims[0],ylims[1])
    
    ## set xlims
    if (xlims is not None):
        plt.xlim(xlims[0],xlims[1])
    
    ## Set savepath
    if (savepath is not None):
        plt.savefig(savepath, dpi=mydpi, bbox_inches='tight')