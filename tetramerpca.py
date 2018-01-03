
__version_info__ = (0, 0, 1)
__version__ = '.'.join(map(str, __version_info__))

import click
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
    'pick':     2,
    'hsp_bit_score_min': 75}
    
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
              help='name of the blast application possibly with path, default is blastn')
            
@click.option('--db',
              default='nt',
              help='blast database name possibly with path, default in nt')
        
@click.option('--num_threads',
              default='12',
              help='blast option passed through, default is 12')
             
@click.option('--num_alignments',
              default='10',
              help='blast option passed through, default is 10')
             
@click.option('--evalue',
              default='10',
              help='blast option passed through, default is 10')    


def main(filename, outdir, window, step, blast_cmd, db, num_threads, 
    num_alignments, evalue, tetra = tetra):
    
    if outdir is None:
        outdir = os.path.dirname(filename)
    
    tetra['filename']           = filename
    tetra['outdir']             = outdir
    tetra['name']               = misc.get_input_name(filename)
    tetra['params']['window']   = window
    tetra['params']['step']     = step
    tetra['params']['hsp_bit_score_min'] = 75
    
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
    
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-window-counts.csv']))
    write_window_counts(tetra['N'], filename = ofile)
    
    # select and write the fails
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-fails.csv']))
    tetra['df_fails'] = tab.tabulate_fails(N, tetra['X'], filename = ofile)
    
    P = PCA(DF)
    tetra['P'] = P
    
    # assemble loadings
    tetra['loadings'] = pd.DataFrame(P.Wt[:, range(0, tetra['params']['npc'])],
        index = list(DF.columns.values), 
        columns =['PC%i' % i for i in range(1,tetra['params']['npc']+1)])    

    # assemble PCs
    tetra['x'] = pd.DataFrame(P.Y[:,range(0, tetra['params']['npc'])], 
        index = list(DF.index.values), 
        columns =['PC%i' % i for i in range(1,tetra['params']['npc']+1)])    

    # select and extract outliers
    tetra['outliers'] = tab.select_outliers_pd(tetra['x'])
    outlier_seqs = misc.extract_outliers(tetra['outliers'], tetra['X'])

    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-counts.csv']))
    ok = inout.write_normalized_counts(tetra['DF'], filename = ofile)
    
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-loadings.csv']))
    ok = inout.write_loadings(tetra['loadings'], filename = ofile)
    
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-tetramer-PC.csv']))
    ok = inout.write_PC(tetra['x'], filename = ofile)
    
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-outliers.fasta']))
    ok = inout.write_fasta(outlier_seqs, filename = ofile)
    
    ifile = ofile
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-outliers.xml']))
    ok = 0 == blast.run_blast(in_file = ifile,  out_file = ofile)
 
    if ok:
        ifile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-outliers.xml']))
        ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-outliers.tsv']))
        nrec = blast.blast_to_table(in_file = ifile, out_file = ofile)
        if nrec > 0:
            tetra['xblast'] = inout.read_blast_table(filename = ofile)   
    
    ofile = os.path.join(tetra['outdir'], ''.join([tetra['name'], '-PC.pdf']))
    r = draw.plot_tetramer(tetra, filename = ofile)
    
    return tetra
    
if __name__=='__main__':
    main()
    