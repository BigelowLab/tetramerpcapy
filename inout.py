import os.path
from Bio import SeqIO
from Bio import SeqUtils
import gzip
import misc
import pandas as pd
import yaml
import io
import pickle

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


def write_params(x, filename = 'example-params.yaml'):
    '''Write a YAML file of the input parameters

    @see https://stackoverflow.com/questions/1773805/how-can-i-parse-a-yaml-file-in-python
    @param x            dict of parameters
    @param filename     fully qualified path
    @return             ummmm
    '''
    with io.open(filename, 'w', encoding='utf8') as outfile:
        yaml.dump(x, outfile, default_flow_style=False, allow_unicode=True)

def read_params(filename = 'example-params.yaml'):
    '''Read a YAML file of the input parameters

    @param filename     fully qualified path
    @return             dict of YAML contents unless the file doesn't exixt then
        default values are returned
    '''

    if os.path.isfile(filename):
        with open(filename, 'r') as stream:
            x = yaml.load(stream)
    else:
        x = {'blastcmd': 'blastn -db nr -num_threads 12 -num_alignments 10 -evalue 10',
            'window':   1600,
            'step':     200,
            'width':    4,
            'npc':      8,
            'pick':     2,
            'hsp_bit_score_min': 75}
    return x


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


def write_fasta(x, filename = 'example-outliers.fasta'):
    '''Write a list of SeqRecord to fasta file

    @param  x a list of SeRecord objects
    @param  filename the file to write to - will overwrite
    @return whatever is returned by Bio.SeqIO.write
    '''
    SeqIO.write(x, filename, "fasta")


def write_window_counts(x, filename = 'example-window-counts.csv'):
    ''' Write a CSV file of window counts per contig.

    @param      dict of contig:count
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''

    x = pd.DataFrame(x.items())
    x.columns = ['cname', 'wcount']
    return x.to_csv(filename)

def read_window_counts(filename = 'example-window-counts.csv'):
    ''' Read a CSV file of window counts per contig.

    @filename   filename and path
    @return     dict contingname:windowcount
    '''

    x = pd.read_table(filename, sep = ',')

    d = dict(zip(x['cname'].values, x['wcount'].values))

    return d



def write_outliers_table(x, filename = 'example-outliers.csv'):
    '''Write a CSV file of selected outliers.

    @param      x a pandas dataframe of selected outliers
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''

    return x.to_csv(filename)

def read_outliers_table(filename = 'example-outliers.csv'):
    '''
    Read a text table of selected outliers.

    @param  filename the name of the file to read
    @return pandas dataframe indexed by contig name
    '''
    x = pd.read_table(filename, sep = ',')
    x = x.set_index('cname', drop = True)

    return x


def write_normalized_counts(x, filename = 'example-counts.csv'):
    '''Write a CSV file of normalized counts per tetramer for each window.

    @param      x a pandas dataframe of normalized counts
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''

    return x.to_csv(filename)

def read_normalized_counts(filename = 'example-counts.csv'):
    '''Read a CSV file of normalized counts per tetramer for each window.

    @filename   filename and path to write to
    @return     pandas dataframe indexed by window name
    '''

    return pd.read_table(filename, sep = ',', index_col = 0)


def write_PC(x, filename = 'example-tetramer-PC.csv'):
    '''Write a CSV file of normalized counts per tetramer for each window.

    @param      x a pandas of PC values
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''
    return x.to_csv(filename)

def read_PC(filename = 'example-tetramer-PC.csv'):
    '''Read a CSV file of normalized counts per tetramer for each window.

    @filename   filename and path to write to
    @return     pandas dataframe indexed by window name
    '''
    return pd.read_table(filename, sep = ',', index_col = 0)

def write_fracs(x, filename = 'example-fracs.pickle'):
    '''Write a PCA fracs dict to file

    @seealos https://pythontips.com/2013/08/02/what-is-pickle-in-python/
    @param x dict
    @param filename the filename and path to write
    @return not sure
    '''
    f = open(filename,'wb')
    pickle.dump(x, f)
    f.close()

def read_fracs(filename = 'example-fracs.pickle'):
    '''Read a PCA fracs dict from file

    @seealos https://pythontips.com/2013/08/02/what-is-pickle-in-python/
    @param filename the filename and path to read
    @return a PCA fracs dict
    '''
    f = open(filename,'rb')
    x = pickle.load(f)
    f.close()
    return x

def write_PCA(x, filename = 'example-PCA.pickle'):
    '''Write a matplotlib.mlab.PCA object to file

    @seealos https://pythontips.com/2013/08/02/what-is-pickle-in-python/
    @param x the matplotlib.mlab.PCA object
    @param filename the filename and path to write
    @return not sure
    '''
    f = open(filename,'wb')
    pickle.dump(x, f)
    f.close()

def read_PCA(filename = 'example-PCA.pickle'):
    '''Read a matplotlib.mlab.PCA object from file

    @seealos https://pythontips.com/2013/08/02/what-is-pickle-in-python/
    @param filename the filename and path to read
    @return matplotlib.mlab.PCA
    '''
    f = open(filename,'rb')
    x = pickle.load(f)
    f.close()
    return x

def write_loadings(x, filename = 'example-tetramer-loading.csv'):
    '''Write a CSV file of loadings per PC for each tetramer.

    @param      x a pandas of PC loading values
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''
    return x.to_csv(filename)


def read_loadings(filename = 'example-tetramer-loading.csv'):
    '''Read a CSV file of normalized counts per tetramer for each window.

    @filename   filename and path to write to
    @return     pandas dataframe indexed by window name
    '''
    return pd.read_table(filename, sep = ',', index_col = 0)

def read_blast_table(filename = 'example-outliers.tsv'):
    '''
    Read a text table of dumped blast XML.

    @param  filename the name of the file to read
    @return pandas dataframe
    '''
    x = pd.read_table(filename, sep = '\t')
    x = x.set_index('Iteration_query_def', drop = False)

    return x


def load_tetramer(path):
    '''Read the elements of TetramerPCA outputs into a dict

    @param path fully qualified path
    @return dict
    '''

    name = os.path.basename(path).split('.')[0]
    tetra = {
        'outdir':   path,
        'filename': os.path.join(path, name + ".fasta.gz"),
        'name':     name}

    # read in either the actual yaml or default values
    tetra['params'] = read_params(filename = name +'-params.yaml')

    # X = input FASTA
    if os.path.exists(tetra['filename']):
        tetra['X'] = read_fasta(filename = tetra['filename'])

    # N = window counts per contig
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-window-counts.csv']))
    if os.path.exists(ifile):
        tetra['N'] = read_window_counts(filename = ifile)

    # df_fails, skip

    # P output of PCA, skip
    #ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-PCA.pickle']))
    #if os.path.exists(ifile):
    #    tetra['P'] = read_PCA(filename = ifile)

    # fracs dict from drawing
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-fracs.pickle']))
    if os.path.exists(ifile):
        tetra['fracs'] = read_fracs(filename = ifile)

    # loadings PCA loadings ('Wt')
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-loadings.csv']))
    if os.path.exists(ifile):
        tetra['loadings'] = read_loadings(filename = ifile)

    # DF normalized counts
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-counts.csv']))
    if os.path.exists(ifile):
        tetra['DF'] = read_normalized_counts(filename = ifile)



    # x PCs
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-tetramer-PC.csv']))
    if os.path.exists(ifile):
        tetra['x'] = read_PC(filename = ifile)


    # xblast
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-outliers.tsv']))
    if os.path.exists(ifile):
        tetra['xblast'] = read_blast_table(filename = ifile)

    # outliers
    ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-outliers.csv']))
    if os.path.exists(ifile):
        tetra['outliers'] = read_outliers_table(filename = ifile)
    return tetra
