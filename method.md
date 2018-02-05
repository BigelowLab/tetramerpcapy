# TetramerPCA Extrema ("Outlier") Selection

### Overview
Tetramers are tallied for one or more contigs.  The frequency of tetramers are
normalized and analyzed using Principle Component Analysis (PCA).  Principle
components are examined pairwise (PC1 v PC2, PC3 v PC4, ...) for extrema.  These
extrema, colloquially outliers, are selected for subsequent analysis using blastn.

### Method
Each contig is divided into windows (typical width 1600 bases).  Overlapping
tetramers are identified as shown below and tallied for that window.  The window
is then advanced one step (typically 200 bases) and the tetramer tally is repeated.

```
AACCAAGGCCACTGGATCCCGAAGGCCTCG
AACC
 ACCA
  CCAA
   CAAG
    AAGG
     ...
```

The table generated is consolidated to fold reverse complements into one tetramer
and then normalized by column as shown below.  Each window is named by contig
name and start-stop window position (1-based) along the contig.

```
cname       wname        AAAA        AAAC        AAAG        AAAT        AACA     ...
    A    A_1-1600 0.011271133 0.011271133 0.006261741 0.006261741 0.008140263     ...
    A  A_201-1800 0.013149656 0.011271133 0.006261741 0.005635567 0.008766437     ...
    A  A_401-2000 0.012523482 0.011271133 0.007514089 0.005009393 0.010644959     ...
    A  A_601-2200 0.012523482 0.011897307 0.008140263 0.005635567 0.011271133     ...
    A  A_801-2400 0.011897307 0.010644959 0.006887915 0.007514089 0.011271133     ...
    A A_1001-2600 0.012523482 0.011897307 0.008140263 0.006887915 0.011897307     ...
    A A_1201-2800 0.010644959 0.009392611 0.009392611 0.006887915 0.010644959     ...
    A A_1401-3000 0.010018785 0.010018785 0.009392611 0.008140263 0.008766437     ...
    A A_1601-3200 0.008766437 0.008766437 0.009392611 0.006261741 0.009392611     ...
    A A_1801-3400 0.006261741 0.007514089 0.008766437 0.006261741 0.008140263     ...
```


Principle components (PC) are computed for the above table; we retain just the
first 8 PCs as shown below (rounded).  Note that we keep record of both the
contig name as well as the window name.

```
cname       wname    PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8
    A    A_1-1600  5.428  0.347  1.267 -0.328  1.273 -0.121 -0.351 -2.254
    A  A_201-1800  5.526  2.187  0.416 -0.486  2.529 -1.207  0.224 -2.084
    A  A_401-2000  3.441  2.361  0.445 -0.147  2.339 -0.810  0.813 -1.983
    A  A_601-2200  3.157  2.179  0.691 -0.295  2.064  0.099  1.291 -2.734
    A  A_801-2400  1.563  2.509  0.603 -0.345  1.236 -0.197  0.382 -3.088
    ...
    I I_1201-2800 -2.680 -0.541 -0.971  1.382  1.118 -0.964 -1.230  0.468
    I I_1401-3000 -1.844 -1.211 -0.984  1.894  1.372 -2.056 -1.129  0.041
    I I_1601-3200 -2.510 -0.719 -0.144  1.437  1.899 -1.934 -0.842 -0.319
    J    J_1-1600 -8.699  3.583 -0.291 -2.268  5.074 -1.796  0.008 -0.145
    J  J_201-1800 -9.249  6.126  0.252 -2.648  5.515 -0.810 -1.026  1.155
```

Next we select extrema (outliers) pairwise (PC1 v PV2, PC3 v PC3, ...) by sorting
the table.  For each pairwise selection, contigs associated with a extreme windows
are removed from subsequent consideration; the table is restored to its initial
state for the next pairwise selection step.  The pseudo-code below shows the process
with examples for a PC1 v PC2 pairwise selection.

```
X = original_PCA_results

for each pairing [PCa and PCb] do

    x = copy(X)  # 11874 rows

    #sort x on PCa which is PC1 in the first iteration
    x = arrange(x, desc(PCa))

    # select most positive windows (PCa's wpos1 and wpos2)
    # cname           wname    PC1   PC2   PC3   PC4    PC5    PC6   PC7   PC8
    #     E     E_4601-6200 13.008 3.999  2.133 2.767 -4.488  0.124 0.277 3.598
    #     A A_650001-651600 11.969 5.379  0.611 2.544 -1.792 -2.043 2.830 2.184


    #remove from x all windows that come from wpos1 and wpos2 contigs
    x = filter(x, cname != "E" & cname != "A")  # leaves 4078 rows

    # select most negative windows (PCa's wneg1 and wneg2)
    #  cname           wname     PC1   PC2   PC3    PC4    PC5    PC6    PC7   PC8
    #      B B_498601-500200 -20.595 1.953 6.705  0.025 -0.934 -0.847 -0.998 2.742
    #      B B_498801-500400 -20.999 2.698 5.572 -0.115 -1.349 -0.679 -0.649 3.150

    #remove from x all windows that come from wneg1 and wneg2 contigs
    x = filter(x, cname != "B") # leaves 1149 rows

    #sort x on PCb which is PC2 in the first iteration
    x = arrange(x, desc(PC2))

    #select most positive windows (PCb's wpos1 and wpos2)
    #     cname       wname   PC1    PC2   PC3    PC4    PC5    PC6    PC7    PC8
    #         F F_4201-5800 3.197 10.869 6.336 -7.581 -0.647 -1.711  0.615  1.301
    #         F F_4001-5600 3.069 10.509 6.133 -7.865 -1.692 -2.556  1.178  1.361

    #remove from x all windows that come from wpos1 and wpos2 contigs
    x = filter(x, cname != 'F') # leaves 1126 rows

    # select most negative windows (PCb's wneg1 and wneg2)
    #  cname           wname    PC1    PC2    PC3    PC4    PC5    PC6   PC7    PC8
    #      C C_143001-144600 -1.402 -5.227 -2.957 -0.776  0.739 -0.206 0.309 -1.361
    #      C   C_94601-96200 -2.943 -5.905  1.164  0.767 -1.242 -0.257 0.525  0.455

```

In this manner we build a list of 32 windows with extreme PCA values. The window
coordinates are then used to subset the original contigs - each window is saved
as its own contig in a FASTA file for use with blastn.
