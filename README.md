For the most up-to-date version of this workflow, please use [tetramers](https://github.com/BigelowLab/tetramers)

## tetramerpcapy

TeteramPCA detects outlier regions in genomes in terms of tetramer composition, which are indicative of either contamination or gene transfer.  Tertamers are
tabulated using an overlapping sliding window. The user specifies both the width of the window and the steps the window is advanced for each tabulation.

This package is adapted from work originally done in perl-language by Alexander Sczyrba at [Bielefeld University Center for Biotechnology](http://www.cebitec.uni-bielefeld.de/~asczyrba) and formerly at [Lawrence Berkeley National Laboratory](http://www.lbl.gov/).

### Requirements:

 In addition to [Python](https://www.python.org) version 3+ and [blast](https://blast.ncbi.nlm.nih.gov) the following packages are required.

 + [Biopython](http://biopython.org/)
 + [click](http://click.pocoo.org/5/)
 + [pandas](http://pandas.pydata.org/index.html)
 + [matplotlib](http://matplotlib.org/1.5.1/index.html)
 + [adjustText](https://pypi.python.org/pypi/adjustText/0.5.3)

### Installation

```
$ cd /path/to/package
$ python setup.py install
```


### Script usage

At a minimum...

```
$ tetramerpca filename
```

... where filename is a fully qualified path to a (possibly gzipped) FASTA file.

```
Usage: tetramerpca [OPTIONS] FILENAME

Options:
  --outdir TEXT          output directory, by default directory same as input
  --window INTEGER       window width in characters, default is 1600
  --step INTEGER         window step size in characters, default is 200
  --blast_cmd TEXT       name of the blast application possibly with path,
                         default is blastn
  --db TEXT              blast database name possibly with path, default in nt
  --num_threads TEXT     blast option passed through, default is 12
  --num_alignments TEXT  blast option passed through, default is 10
  --evalue TEXT          blast option passed through, default is 10
  -h, --help             Show this message and exit.

```

### What the process does

 + Generates a matrix of tetramer frequencies (columns) of each contig window (rows), using the window and step sizes. The number of tetramers is 136, reduced from the maximum of 256 to eliminate reverse-complementary tetramers.

 + Normalizes the tetramer frequency table as an input to estimate the first eight principal components.

 + For each of the principal components, identifies those contigs that contain "outliers", i.e. windows with either extreme positive or negative values.

 + Generate blastn output for the selected outliers.

 + Generates four scatter-plots: PC1 vs PC2, PC3 vs PC4, PC5 vs PC6 and PC7 vs PC8. All dots are grey except for those corresponding to contigs that contain outlier windows. Each outlier-containing contig is plotted using separate color and/or symbol, so that they are clearly visible on the plot.  A color/symbol legend is generated for each plot.

 + Output a variety of files where `name` is from the input file.

   - Normalized tetramer frequency table, `name-counts.csv`

   - Principal components and loadings, `name-loading.csv`

   - PCA plots, `name-outliers.pdf`

   - Outlier selection table, `name-outliers.csv`

   - Outlier window contigs, `name-outliers.fasta`

   - blastn output, `name-outliers.xml`

   - blastn output as a table, `name-outliers.tsv`


   - a listing of failed contigs, `name-failed.csv`

### Loading existing output

For the programmer's convenience existing output from a prior run can be loaded into a dictionary.

```
import tetramerpca as tet
tetra = tet.inout.load_tetramer('/path/to/the/exisiting/output/data/set')
```
