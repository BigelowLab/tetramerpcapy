import os.path
from Bio import SeqIO
from Bio import SeqUtils
import gzip
import misc
import pandas as pd

def is_gzipped(filename = 'example.fasta.gz'):

    '''Test input file to determine if it is gzipped.
    
    @param filename     fully qualified path
    @return boolean     true if the file is gzipped
    '''
    try:
        file_handle = open(filename,"rb")
    except IOerror, e:
        print 'error opening file:', e
        
    x1 = file_handle.read(1)
    x2 = file_handle.read(1)
    file_handle.close()

    is_gz = (x1 == '\x1f') and (x2 == '\x8b')
    
    return is_gz
    
    
def read_fasta(filename = 'example.fasta.gz'):
    '''Read a fasta file (possibly gzipped)
    
    @param filename     fully qualified path
    @return             list of zero or more SeqRecord
    '''
    contigs = list()
    
    if not os.path.isfile(filename):
        print 'file not found:', filename
        return contigs
    
    if is_gzipped(filename):
        handle = gzip.open(filename, "r")
    else :
        handle = open(filename, "r")
        
    contigs=list(SeqIO.parse(handle, 'fasta'))
    
    handle.close()
    
    return contigs
    
    
def write_normalized_counts(x, filename = 'example-counts.csv'):
    '''Write a CSV file of normalized counts per tetramer for each window.
    
    @param      x a pandas dataframe of normalized counts
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''
    
    return x.to_csv(filename)
    
def write_PC(x, filename = 'example-tetramer-PC.csv'):
    '''Write a CSV file of normalized counts per tetramer for each window.
    
    @param      x a pandas of PC values
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''    
    return x.to_csv(filename)

def write_loadings(x, filename = 'example-tetramer-loading.csv'):
    '''Write a CSV file of loadings per PC for each tetramer.
    
    @param      x a pandas of PC loading values
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''    
    return x.to_csv(filename)
    
def write_fasta(x, filename = 'example-outliers.fasta'):
    '''Write a list of SeqRecord to fasta file
    
    @param  x a list of SeRecord objects
    @param  filename the file to write to - will overwrite
    @return whatever is returned by Bio.SeqIO.write
    '''
    SeqIO.write(x, "example-outliers.fasta", "fasta")
    
    
def read_blast_table(filename = 'example-outliers.tsv'):
    '''
    Read a text table of dumped blast XML.

    @param  filename the name of the file to read
    @return pandas dataframe
    '''
    x = pd.read_table(filename, sep = '\t')
    x = x.set_index('Iteration_query_def', drop = False)

    return x