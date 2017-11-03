
__version_info__ = (0, 0, 1)
__version__ = '.'.join(map(str, __version_info__))

import Bio
import numpy as np
import pandas as pd
from matplotlib.mlab import PCA
import os

import inout
import lut
import misc
import tabulate as tab
import blast
import draw

tetra = dict()

tetra['params'] = {
    'blastcmd': 'blastn -db nr -num_threads 12 -num_alignments 10 -evalue 10',
    'window':   1600,
    'step':     200,
    'width':    4,
    'npc':      8,
    'pick':     2}
    
@click.command(context_settings=dict(help_option_names=['-h', '--help']))

@click.argument('filename')

@click.option('--outdir',
              default=None,
              help='output directory, by default directory same as input')

@click.option('--window',
              default=tetra['params']['window'],
              type=click.INT,
              help='window width in characters, default is 1600')

@click.option('--step',
              default=tetra['params']['step'],
              type=click.INT,
              help='window step size in characters, default is 200')

@click.option('--blast_cmd',
              default='blastn',
              help='name of the bast application possibly with path')
            
@click.option('--db',
              default='nt',
              help='blast database name possibly with path')
        
@click.option('--num_threads',
              default='12',
              help='blast option passed through')
             
@click.option('--num_alignments',
              default='10',
              help='blast option passed through')
             
@click.option('--evalue',
              default='10',
              help='blast option passed through')    


def main(filename, outdir, window, step, blast_cmd, db, num_threads, 
    num_alignments, evalue, tetra):
    
    if outdir is None:
        outdir = os.path.dirname(filename)
    
    tetra['filename']           = filename
    tetra['name']               = misc.get_input_name(filename)
    tetra['params']['window']   = window
    tetra['params']['step']     = step
    
    fmt = "%s -dd %s -num_threads %s -num_alignments %s -evalue %s"
    tetra['params']['blastcmd']=  fmt % (blast_cmd, db, num_threads, 
        num_alignments, evalue)
    
    tetra['X'] = inout.read_fasta(filename = filename)
    if len(tetra['X']) == 0:
        return tetra
        
    # tally the tetramers
    N, DF = tab.tabulate_tetramers(tetra['X'])
    tetra['N'] = N
    tetra['DF'] = DF
    
    # select and write the shorties
    tetra['df_fails'] = tab.tabulate_fails(N, tetra['X'])
    
    # obvious
    P = PCA(DF)
    tetra['P'] = P
    
    # assemble and write the loadings
    tetra['loadings'] = pd.DataFrame(P.Wt[:, range(0, tetra['params']['npc'])],
        index = list(DF.columns.values), 
        columns =['PC%i' % i for i in range(1,tetra['params']['npc']+1)])    

    # assemble and write the PCs
    tetra['x'] = pd.DataFrame(P.Y[:,range(0, tetra['params']['npc'])], 
        index = list(DF.index.values), 
        columns =['PC%i' % i for i in range(1,tetra['params']['npc']+1)])    

    # select and extract outliers
    tetra['outliers'] = tab.select_outliers(tetra['x'])
    outlier_seqs = misc.extract_outliers(tetra['outliers'], tetra['X'])

    ok = inout.write_normalized_counts(tetra['DF'], 
        filename = 'example-counts.csv')
    ok = inout.write_loadings(tetra['loadings'], 
        filename = 'example-loadings.csv')
    ok = inout.write_PC(tetra['x'], 
        filename = 'example-tetramer-PC.csv')
    ok = inout.write_fasta(outlier_seqs,
        filename = 'example-outliers.fasta')
    
    ok = 0 == blast.run_blast(in_file = 'example-outliers.fasta', 
        out_file = 'example-outliers.xml')
 
    if ok:
         nrec = blast.blast_to_table(in_file = 'example-outliers.xml',
             out_file = 'example-outliers.tsv')
         if nrec > 0:
             tetra['xblast'] = inout.read_blast_table(filename = 'example-outliers.tsv')   
    
    r = draw.plot_tetramer(tetra, filename = 'example-PC.pdf')
    
    return tetra
    
if __name__=='__main__':
main()
    