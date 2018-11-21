import pandas as pd
import string
import os.path

def get_input_name(filename = 'example.fasta.gz'):
    '''
    Return the extensionless file basename.  Assumes that any . indicates extension.
    
    foo         yields foo
    foo.bar     yields foo
    foo.bar.gz  yields foo
    
    @param filename the name of the file
    @return the extensionless basename
    '''
    bname = os.path.basename(filename).split('.')[0]
    return bname

def create_empty_dataframe(columns = [('name', str), ('score', int),
        ('height', float), ('weight', float)]):
    ''' Create an empty DataFrame with the specified columns
    
    See https://stackoverflow.com/questions/38523965/python-pandas-create-empty-dataframe-specifying-column-dtypes

    @param columns list of column name-type tuples
    @return an empty DataFrame
    '''
    
    index = pd.Index([], name="id", dtype=int)
    # create the dataframe from a dict
    return pd.DataFrame({k: pd.Series(dtype=t) for k, t in columns})


def split_name(nm = 'foo_a_b_1-500'):
    ''' Split the input into 4 parts 
    
    [name (foo_a_b), window (1-500), start (1) end (500)]
    
    @param nm the string to split
    @return a four element list 
    '''
    # a = string.rsplit(nm, "_", 1)
    # b = string.split(a[1], "-")
    a = nm.rsplit("_", 1)
    b = a[1].split("-")
    return a + [int(v) - 1 for v in b]  

def parse_windowname(wn = ['foo_a_b_1-500', 'foo_a_b_501-1000']):
    '''Parse a list of window names into conting name, window id
    
    foo_bah_baz-1-500 is split into foo_bar_baz and 1-500
    
    @param  wn list of window names
    @return pandas DataFrame
    '''
    if not isinstance(wn, list):
        if isinstance(wn, str):
            wn = [wn]
        else:
            wn = list(wn)
           
    s = [split_name(nm = s) for s in wn]
    
    df = pd.DataFrame(s,
        columns = ['name', 'window', 'start', 'stop'],
        index = wn)
    
    return df
    
 
def extract_seq(s, nv):
    '''Extract a subsequence from a sequence record
    
    @param s a SeqRecord
    @param nv pandas Series with name, window, start and stop
    @return slices SeqRecord
    '''
    
    ss = s[int(nv.start):(int(nv.stop)+1)]
    ss.id = '{0}_{1}'.format(nv['name'], nv['window'])
    # copy the above so they are not written to the header line
    ss.name = ss.id
    ss.description = ss.id
    
    return ss
           
def extract_outliers(x, X):
    '''Create a list of SeqRecords
    
    @param      x a pandas of outlier table
    @param      X list of one or more Bio.SeqRecord
    @filename   filename and path to write to
    @return     whatever is returned by pandas dataframe.to_csv()
    '''
 
    # convert the SeqRecord list to a dict 
    XID = [s.id for s in X]
    X = dict(zip(XID, X))
    
    ss = list()
    #for contig in list(x.index.values):
    #    #one or more contigs may be involved
    #    wn = x.loc[contig,'name']
    #    #if isinstance(wn, str):
    #    #    wn = [wn]
    #    #else:
    #    #    wn = list(wn)
    #    nav = parse_windowname(wn)
    #    ss = ss + [extract_seq(X[contig], row) for index, row in nav.iterrows()] 
    
    # the above doesn't iterate over rows, but over row names which is muddles
    # the loops, instead we iterate using iterrows
    # with returns an iteration count and a <something> whic contains the data
    # in the row
    # because of the way I established extract_seq it expects a row of a
    # dataframe so we user iterrows a second time for a 1-row dataframe.  
    # Could be better.
    for i,r in x.iterrows():
        nav = parse_windowname(r['name'])
        ss = ss + [extract_seq(X[nav['name'][0]], row) for index, row in nav.iterrows()]
        
    return ss


