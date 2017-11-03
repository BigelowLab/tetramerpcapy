# for each seq
#   tabulate tetramers (and comps) in a Window
#   advance the Window by Step and retabulate

import lut
import numpy as np
import pandas as pd
import misc
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqUtils
from collections import OrderedDict


def get_range(y, col = 'PC1'):
    '''Given a indexed dataframe find the max/min windows and values.

    @param y    pandas dataframe with rownames
    @param col  the column to operate upon
    @return     dataframe indexed by contig name
    '''

    #indices low to high

    imin = y.loc[:,col].idxmin()
    imax = y.loc[:,col].idxmax()
    r = pd.DataFrame({'min_name':  imin, 
                      'min_value': y.loc[imin,col], 
                      'max_name':  imax, 
                      'max_value': y.loc[imax,col]},
                      index = [imin.rsplit("_",1)[0]] )
    return r

def select_outliers(x, pick = 2):
    '''Select the outliers from the PCA results. 

    This is a complex arrangement of simple steps to identify the 'extreme' 
    response windows within a set of contigs in each PC direction. Extreme means
    here most positive ('hi') and most negative ('lo') windows. Selection 
    is done in pairs of PC's (PC1 + PC2, PC3 + PC4, ...).  Within each pairing
    the selection set is winnowed to prevent the selection of a given contig
    more than once - the idea being to expand the pool of blast candidates. 

    @param x    a dataframe of PC1, PC2, ..., PC8 indexed by window name
    @param pick the number of extreme candidates to select.
    @return     dataframe
    '''

    nmpc = x.columns.tolist()
    nm = misc.parse_windowname(list(x.index.values))

    xx = dict(list(x.groupby(nm['name'])))
    pick_hi_idx = range(0, pick)
    pick_lo_idx = pick_hi_idx[::-1]
    hi_names = ['hi'] * pick
    lo_names = ['lo'] * pick
    R = list()
    for i in range(0, len(nmpc),2):
        pcname = nmpc[i]
        keys = xx.keys()
        # this provides an contigname-indexed dataframe
        r = pd.concat([get_range(xx[k], col = nmpc[i]) for k in keys])

        # select the two most positive hits
        r = r.sort_values('max_value', ascending = False)
        hi1_contig = r.iloc[pick_hi_idx][['max_name']]
        hi1_contig['typ'] = hi_names
        hi1_contig['PC']  = [pcname] * pick
        hi1_contig.columns = ['name', 'typ', 'PC']

        # now exclude those selected above and 
        # select the two most negative hits
        r = r.drop(list(hi1_contig.index))
        r = r.sort_values('min_value', ascending = True)
        lo1_contig = r.iloc[pick_lo_idx][['min_name']]
        lo1_contig['typ'] = lo_names
        lo1_contig['PC']  = [pcname] * pick
        lo1_contig.columns = ['name', 'typ', 'PC']

        R1 = pd.concat([hi1_contig, lo1_contig])

        # now we do the next PC - note we do this is pairs so that we keep 
        # trimming the selection between PC1 -> PC2, then we start with the
        # fullset for PC3 -> PC4.  I can't recall why we do this but it must have
        # been important.

        # remove the keys we have already used
        remove_contigs = list(R1.index)
        dummy = [keys.remove(j) for j in remove_contigs]

        pcname = nmpc[i+1]

        r = pd.concat([get_range(xx[k], col = nmpc[i+1]) for k in keys])
        # select the two most positive hits
        r = r.sort_values('max_value', ascending = False)
        hi2_contig = r.iloc[pick_hi_idx][['max_name']]
        hi2_contig['typ'] = hi_names
        hi2_contig['PC']  = [pcname] * pick
        hi2_contig.columns = ['name', 'typ', 'PC']

        # now exclude those selected above and 
        # select the two most negative hits
        r = r.drop(list(hi2_contig.index))
        r = r.sort_values('min_value', ascending = True)
        lo2_contig = r.iloc[pick_lo_idx][['min_name']]
        lo2_contig['typ'] = lo_names
        lo2_contig['PC']  = [pcname] * pick
        lo2_contig.columns = ['name', 'typ', 'PC']

        R2 = pd.concat([hi2_contig, lo2_contig])        

        R.append(pd.concat([R1,R2]))

    R = pd.concat(R)
    return R



def tabulate_fails(N, X, filename = None):

    '''Tabulates failed contigs (no windows) and write to CSV

    @param N    dict of window counts by contig
    @param X    the SeqRecords database (fasta file origin)
    @param      filename the file to write to or Noen to ignore
    @return     pandas DataFrame of failed contigs
    '''

    z = []
    nm = []
    status = 'short'
    for x in X:
        n = N[x.id]
        if n <= 0:
            z.append([ status,len(x), SeqUtils.GC(x.seq) / 100])
            nm.append(x.id)

    df = pd.DataFrame(z, 
        columns = ['status', 'contigLength', 'pGC'], 
        index = nm)

    if filename is not None:
        df.to_csv(filename)

    return df


def get_dictnames(dct = lut.get_utetramers(), sort_it = 'ascending'):
    '''Retrieve names of dictionary, possibly sorted

    @param x        the dictionary
    @param sort_it  string either ascending or descending
    @return         list of possibly sorted keys
    '''

    k = dct.keys()
    if sort_it == 'ascending':
        k.sort()
    elif sort_it == 'descending':
        k.sort(reverse = True)

    return k



def tabulate_kmers(v = ['GGTT', 'GTTA', 'TTAC', 'TACA', 'ACAN', 'CANN', 'ANNN'], 
        kmers = get_dictnames()):
    '''Tabulate the kmers in v by the elements in kmers

    @params v       list of items that will be counted
    @parmas kmers   list of items to count 
    @return         list of counts - in the order of kmers
    '''
    x = [v.count(kmers[i]) for i in xrange(len(kmers)) ] 

    return x


def extract_kmers(x = 'AATTTTTTTCGTTCTGCGGAACCACCAAAACCAGTC', k = 4 , step = 1):
    '''Chop a string into k-mers advancing by step.

    If step < k then this performs a rolling chop with characters being repeated
    between successive kmers.

    @param x    the input string
    @param k    the width of kmers in charcaters
    @param step the number of characters to advance for each kmer
    @return     list of kmers
    '''

    kmers = [x[i:i+k] for i in xrange(0, len(x) - k + 1, step)]

    return kmers


def tabulate_seq(x, window = 1600, step = 200, width = 4):
    '''Tabulate tetramers in a Bio.Seq

    @param window   the width of the sliding window in characters
    @param step     the number of characters to advance the window each iteration
    @param width    the width of the kmers to tally
    @return         a pandas.dataframe indexed by window name ala id-start-stop or None
    '''

    DNALETTERS = lut.get_dnaletters()
    UTETRAMERS = lut.get_utetramers()
    tetnm = UTETRAMERS.keys()
    tetnm.sort()

    # sequential start locations - note that we want a np.ndarray
    # because we will increment the values by adding a scaler which is not
    # as efficient when done with a list (which requires explicit iteration)
    starts = np.arange(window - width)
    nstarts = len(starts)
    nx = len(x)

    if nx >= window:
        nsteps = ( ( nx - window ) / step ) + 1
        df = pd.DataFrame(0, index = np.arange(nsteps), columns = tetnm)
        rownames = [' ']*nsteps

        for istep in np.arange(nsteps):
            s0 = istep * step
            s1 = s0 + starts
            rownames[istep] = '{0}_{1}-{2}'.format(x.id, s1[0] + 1, s1[nstarts-1] + width + 1)
            s = str(x[s1[0]:s1[nstarts-1] + width + 1].seq)
            df.ix[istep] = tabulate_kmers(extract_kmers(s), kmers = tetnm)

        df.index = rownames 
    else:
        df = None;

    return df

def tabulate_seqs(X, window = 1600, step = 200, width = 4):
    '''Tabulate tetramers in a Bio.SeqRecord

    @param window   the width of the sliding window in characters
    @param step     the number of characters to advance the window each iteration
    @param width    the width of the kmers to tally
    @return         a list of pandas.dataframe as produced by tabulate_seq
    '''

    # iterate through them all
    dfs = [tabulate_seq(x, window = window, step = step, width = width) for x in X]

    return dfs

def reduce_tetramers(x):
    
    comp = lut.get_rcomp_tetramers()
    keep = lut.get_keep_tetramers()

    # first - to the keep columns add the rcomp columns
    for k in comp.keys():
        x[k] = x[k] + x[comp[k]]

    # next subset the columns, retaining only the keepers
    x = x.loc[:,keep]

    return x

def normalize_tetramers(x):
    ''' Transform a dataframe of counts to normalized (0-1) counts.

    Rows that are everywhere 0 (i.e. no teteramers found) are set to zero everywhere
    
    @param x    pandas dataframe of tetramer counts
    @return     pandas dataframe of same shape with normalize counts
    '''

    s = x.sum(axis = 1)
    # see https://stackoverflow.com/questions/26537878/pandas-sum-across-columns-and-divide-each-cell-from-that-value
    x = (x.T / x.T.sum()).T
    x = x.fillna(0)

    return x

def tabulate_tetramers(X, window = 1600, step = 200, width = 4):
    '''Tabulate contigs in a Bio.SeqRecord

    Returns a tuple of 
        (1) counts of windows by contig (even failed ones which will be 0)
        (2) dataframe of tetramer counts by tetramer (cols) and windows (rows)

    @param X        
    @param window   the width of the sliding window in characters
    @param step     the number of characters to advance the window each iteration
    @param width    the width of the kmers to tally
    @return         a tuple
    '''

    dfs = tabulate_seqs(X, window = window, step = step, width = width)

    # create a counts dictionary - I miss R
    id = [x.id for x in X]
    n = [0 if (x is None) else x.shape[0] for x in dfs]
    counts = dict(zip(id, n))

    return counts, normalize_tetramers(reduce_tetramers(pd.concat(dfs)))
