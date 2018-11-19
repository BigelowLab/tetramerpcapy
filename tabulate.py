# for each seq
#   tabulate tetramers (and comps) in a Window
#   advance the Window by Step and retabulate

import lut
import numpy as np
import pandas as pd
import misc
import string
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqUtils
from collections import OrderedDict
from collections import Counter


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


def select_outliers_pd_permissive(x, unpack = True):
    '''
    Selects 'outliers' from a [nwindow x nPC] dataframe using permissive reductions
    
    index the input with cn (contigname) and wn (windowname)

    for each PC do:

      # extract the column PC and sort
      # select the two most positive (wpos)
      # filter out matching contigs
      # if remaining is < 2 rows then add the wpos contigs back in
      # select the two most negative (wneg)

      pc = x.loc[:, PC].copy().sort_values(inplace = True)
      idx = m.index.values
      n = len(idx)
      wpos = [x[1] for x in idx[(n-2):]]

      

    @param x dataframe [nwindows x nPC], indexed by windowname (wname)
    @param unpack boolean, if True reshape to match the tetramerpca workflow
    @return dataframe depending upon unpack value
      unpack = False, indexed by PC with 4 columns wpos1,wpos2, wneg1, wneg2 where
        the columns are the window names ala "contigname_start-stop"
      unpack = True, indexed by contigname (cname) with 3 columns typ (hi/lo),
        PC (PC1, PC2, ...) and name (wname)
    '''
    def unpack_select(x):
        '''
        Given the dense outlier selection data frame, make the long form used by
        tetramerpca.

        @param x a data frame indexed by PC wih 4 columns wpos1, wpos2, wneg1 and wneg2
        @return a data frame indexed by contig name with columns name (wname), typ,
            and PC
        '''
        x = x.unstack().reset_index()
        x = x.replace("wpos1", "hi")
        x = x.replace("wpos2", "hi")
        x = x.replace("wneg1", "lo")
        x = x.replace("wneg2", "lo")
        x.columns = ['typ', 'PC', 'name']
        x['cname'] = extract_cname(list(x['name'].values))
        x = x.set_index('cname')
        return x

    def split_name(nm = 'foo_a_b_1-500'):
        ''' Split the input into 4 parts

        [name (foo_a_b), window (1-500), start (1) end (500)]

        @param nm the string to split
        @return a four element list
        '''
        #a = string.rsplit(nm, "_", 1)
        a = nm.rsplit("_", 1)
        #b = string.split(a[1], "-")
        b = a[1].split("-")
        return a + [int(v) - 1 for v in b]

    def extract_cname(wn = ["A_1-1600","B_201-1800"]):
        if not isinstance(wn, list):
            if isinstance(wn, str):
                wn = [wn]
            else:
                wn = list(wn)
        s = [split_name(nm = s)[0] for s in wn]
        return s


    wname = list(x.index.values)
    cname = extract_cname(wname)
    x.set_index([cname, wname], inplace = True)
    x.index.set_names(['cn', 'wn'], inplace = TRUE)

    r = pd.DataFrame(data = 'foo',
        index = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8'],
        columns=['wpos1','wpos2','wneg1','wneg2'])

    for PC in x.columns.values.tolist():
        pc = x.loc[:, PC].copy()
        pc.sort_values(inplace = True)
        ipos = pc.tail(2).index.values
        wpos = [x[1] for x in ipos]
        cpos = [x[0] for x in ipos]
        r.at[PC, 'wpos1'] = wpos[1]
        r.at[PC, 'wpos2'] = wpos[0]
        n = Counter(~pc.index.isin(cpos, level = 'cn'))
        if n[True] >= 2:
            pc.drop(cpos, level = 'cn', inplace = True)
        ineg = pc.head(2).index.values
        wneg = [x[1] for x in ineg]
        r.at[PC, 'wneg1'] = wneg[0]
        r.at[PC, 'wneg2'] = wneg[1]
        
    x.reset_index(inplace = True, drop = True)
    x.set_index(wname, inplace = True)
    
    if unpack:
        r = unpack_select(r)

    return r


def select_outliers_pd(x, unpack = True, permissive = True):
    '''
    Selects 'outliers' from a [nwindow x nPC] dataframe.
    Selection of outliers occurs in PCa vs PCb (e.g. PC1 v PC2, PC3 v PC4, ...)
    pairs.  Here is pseudocode that operates on input_df dataframe

    for each PCa vs PCb pairing do:

      df = input_df.copy()

      sort df on PCa
      select most positive windows (PCa's wpos1 and wpos2)
      remove from df all windows that come from wpos1 and wpos2 contigs
      select most negative windows (PCa's wneg1 and wneg2)
      remove from df all windows that come from wneg1 and wneg2 contigs

      sort df on PCb
      select most positive windows (PCb's wpos1 and wpos2)
      remove from df all windows that come from wpos1 and wpos2 contigs
      select most negative windows (PCb's wneg1 and wneg2)


    @param x dataframe [nwindows x nPC], indexed by windowname (wname)
    @param unpack boolean, if True reshape to match the tetramerpca workflow
    @param permissive boolean, if True then use permissive filtering
    @return dataframe depending upon unpack value
      unpack = False, indexed by PC with 4 columns wpos1,wpos2, wneg1, wneg2 where
        the columns are the window names ala "contigname_start-stop"
      unpack = True, indexed by contigname (cname) with 3 columns typ (hi/lo),
        PC (PC1, PC2, ...) and name (wname)
    '''

    if permissive:
        return(select_outliers_pd_permissive(x, unpack = unpack))


    def unpack_select(x):
        '''
        Given the dense outlier selection data frame, make the long form used by
        tetramerpca.

        @param x a data frame indexed by PC wih 4 columns wpos1, wpos2, wneg1 and wneg2
        @return a data frame indexed by contig name with columns name (wname), typ,
            and PC
        '''
        x = x.unstack().reset_index()
        x = x.replace("wpos1", "hi")
        x = x.replace("wpos2", "hi")
        x = x.replace("wneg1", "lo")
        x = x.replace("wneg2", "lo")
        x.columns = ['typ', 'PC', 'name']
        x['cname'] = extract_cname(list(x['name'].values))
        x = x.set_index('cname')
        return x

    def split_name(nm = 'foo_a_b_1-500'):
        ''' Split the input into 4 parts

        [name (foo_a_b), window (1-500), start (1) end (500)]

        @param nm the string to split
        @return a four element list
        '''
        #a = string.rsplit(nm, "_", 1)
        a = nm.rsplit("_", 1)
        #b = string.split(a[1], "-")
        b = a[1].split("-")
        return a + [int(v) - 1 for v in b]

    def extract_cname(wn = ["A_1-1600","B_201-1800"]):
        if not isinstance(wn, list):
            if isinstance(wn, str):
                wn = [wn]
            else:
                wn = list(wn)
        s = [split_name(nm = s)[0] for s in wn]
        return s

    orig_x = x.copy()
    wname = list(orig_x.index.values)
    orig_x.reset_index(0, drop = True, inplace = True)
    orig_x.insert(0, 'wname', wname)
    cname = extract_cname(wname) # extract_cname(list(x['wname'].values))
    orig_x.insert(0, 'cname', cname)
    orig_x = orig_x.set_index("cname")
    r = pd.DataFrame(data = 'foo',
        index = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8'],
        columns=['wpos1','wpos2','wneg1','wneg2'])

    for iPC in list(range(1,8,2)):

        # start the paired selection with a fresh copy of input PC matrix
        # sorted by the PCa values
        PCa = "%s%i" % ('PC', iPC)
        PCb = '%s%i' % ('PC', iPC + 1)

        x       = orig_x.copy().sort_values(by = PCa)
        n       = x.shape[0]
        whi     = x.iloc[n-1, 0]                    # first place positive
        chi     = extract_cname(whi)[0]
        x       = x.drop(chi)                      # remove contigs
        n       = x.shape[0]
        whi     = [whi, x.iloc[n-1, 0]]             # [first,second] place positive
        chi     = extract_cname(whi)
        x       = x.drop(chi)                       # remove contigs
        wlo     = x.iloc[0, 0]                      # first place negative
        clo     = extract_cname(wlo)[0]
        x       = x.drop(clo)                       # remove contigs
        wlo     = [wlo, x.iloc[0, 0]]               # [first, second] place negative

        # save [most pos, pos, neg, most neg]
        r.loc[PCa] = [whi[0], whi[1], wlo[1], wlo[0]]

        x       = x.drop(extract_cname(wlo))        # remove contigs
        x       = x.sort_values(by = PCb)           # second PC
        n       = x.shape[0]
        whi     = x.iloc[n-1, 0]                    # first place positive
        chi     = extract_cname(whi)
        x       = x.drop(chi)                       # remove contigs
        n       = x.shape[0]
        whi     = [whi, x.iloc[n-1, 0]]             # [first,second] place positive
        chi     = extract_cname(whi)
        x       = x.drop(chi)                       # remove contigs
        wlo     = x.iloc[0, 0]
        clo     = extract_cname(wlo)[0]
        x       = x.drop(clo)                       # remove contigs
        wlo     = [wlo, x.iloc[0, 0]]               # [first, second] place negative

        # save those  [most pos, pos, neg, most neg]
        r.loc[PCb] = [whi[0], whi[1], wlo[1], wlo[0]]

    if unpack:
        r = unpack_select(r)

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
    pick_hi_idx = list(range(0, pick))
    pick_lo_idx = pick_hi_idx[::-1]
    hi_names = ['hi'] * pick
    lo_names = ['lo'] * pick
    R = list()
    for i in list(range(0, len(nmpc),2)):
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

    # k = dct.keys()
    # https://stackoverflow.com/questions/22845330/attribute-error-in-python-wont-go-away#22845343
    k = list(dct.keys())
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
    x = [v.count(kmers[i]) for i in list(range(len(kmers))) ]

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

    kmers = [x[i:i+k] for i in list(range(0, len(x) - k + 1, step))]

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
    tetnm = list(UTETRAMERS.keys())
    tetnm.sort()

    # sequential start locations - note that we want a np.ndarray
    # because we will increment the values by adding a scaler which is not
    # as efficient when done with a list (which requires explicit iteration)
    starts = np.arange(window - width)
    nstarts = len(starts)
    nx = len(x)

    if nx >= window:
        nsteps = int(( ( nx - window ) / step ) + 1)
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
