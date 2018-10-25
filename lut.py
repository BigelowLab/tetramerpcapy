
def get_dnaletters():
    ''' Return a dictionary of DNA letters
    
    @return dictionary
    '''
    return {'A':1, 'C':2, 'G':3, 'T':4}
    
 
def get_rcomp_tetramers():
    ''' Retrieve a 120 element dict of reverse complement teteramers
    
    @return dict of tetramer= reverse_comp pairs
    '''
    
    tet = ["AAGT", "ACCT", "AACT", "AGAT", "ACAT", "AAAT", "ATTG", "AGTG", 
        "ACTG", "AATG", "ATGG", "AGGG", "ACGG", "CAGG", "AAGG", "ATCG", 
        "AGCG", "CCCG", "ACCG", "CACG", "AACG", "ATAG", "CGAG", "AGAG", 
        "CCAG", "ACAG", "CAAG", "AAAG", "CTTC", "ATTC", "CGTC", "AGTC", 
        "CCTC", "ACTC", "CATC", "AATC", "CTGC", "ATGC", "CGGC", "AGGC", 
        "CCGC", "ACGC", "GAGC", "CAGC", "AAGC", "CTCC", "ATCC", "CGCC", 
        "AGCC", "GCCC", "CCCC", "ACCC", "GACC", "CACC", "AACC", "CTAC", 
        "ATAC", "GGAC", "CGAC", "AGAC", "GCAC", "CCAC", "ACAC", "GAAC", 
        "CAAC", "AAAC", "GTTA", "CTTA", "ATTA", "GGTA", "CGTA", "AGTA", 
        "GCTA", "CCTA", "ACTA", "GATA", "CATA", "AATA", "GTGA", "CTGA", 
        "ATGA", "GGGA", "CGGA", "AGGA", "GCGA", "CCGA", "ACGA", "TAGA", 
        "GAGA", "CAGA", "AAGA", "GTCA", "CTCA", "ATCA", "GGCA", "CGCA", 
        "AGCA", "TCCA", "GCCA", "CCCA", "ACCA", "TACA", "GACA", "CACA", 
        "AACA", "GTAA", "CTAA", "ATAA", "TGAA", "GGAA", "CGAA", "AGAA", 
        "TCAA", "GCAA", "CCAA", "ACAA", "TAAA", "GAAA", "CAAA", "AAAA"]
    
    return rcomp_kmers(tet)
    
def get_keep_tetramers():
    ''' Retrieve a list of tetramers to retain
    
    @param a list of 136 tetramers
    '''
    
    keepers = [
        "AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", 
        "AAGA", "AAGC", "AAGG", "AAGT", "AATA", "AATC", "AATG", "AATT", 
        "ACAA", "ACAC", "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", 
        "ACGA", "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "AGAA", 
        "AGAC", "AGAG", "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", "AGGA", 
        "AGGC", "AGGG", "AGTA", "AGTC", "AGTG", "ATAA", "ATAC", "ATAG", 
        "ATAT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATGG", "ATTA", 
        "ATTC", "ATTG", "CAAA", "CAAC", "CAAG", "CACA", "CACC", "CACG", 
        "CAGA", "CAGC", "CAGG", "CATA", "CATC", "CATG", "CCAA", "CCAC", 
        "CCAG", "CCCA", "CCCC", "CCCG", "CCGA", "CCGC", "CCGG", "CCTA", 
        "CCTC", "CGAA", "CGAC", "CGAG", "CGCA", "CGCC", "CGCG", "CGGA", 
        "CGGC", "CGTA", "CGTC", "CTAA", "CTAC", "CTAG", "CTCA", "CTCC", 
        "CTGA", "CTGC", "CTTA", "CTTC", "GAAA", "GAAC", "GACA", "GACC", 
        "GAGA", "GAGC", "GATA", "GATC", "GCAA", "GCAC", "GCCA", "GCCC", 
        "GCGA", "GCGC", "GCTA", "GGAA", "GGAC", "GGCA", "GGCC", "GGGA", 
        "GGTA", "GTAA", "GTAC", "GTCA", "GTGA", "GTTA", "TAAA", "TACA", 
        "TAGA", "TATA", "TCAA", "TCCA", "TCGA", "TGAA", "TGCA", "TTAA" ]
    
    return keepers
    
     
def get_utetramers():
    ''' Retrieve a dictionary of 256 tetramer combinations
    
    @return dictionary
    '''   
    utetramers = [
        "AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", 
        "AAGA", "AAGC", "AAGG", "AAGT", "AATA", "AATC", "AATG", "AATT", 
        "ACAA", "ACAC", "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", 
        "ACGA", "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", 
        "AGAA", "AGAC", "AGAG", "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", 
        "AGGA", "AGGC", "AGGG", "AGGT", "AGTA", "AGTC", "AGTG", "AGTT", 
        "ATAA", "ATAC", "ATAG", "ATAT", "ATCA", "ATCC", "ATCG", "ATCT", 
        "ATGA", "ATGC", "ATGG", "ATGT", "ATTA", "ATTC", "ATTG", "ATTT", 
        "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT", 
        "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT", 
        "CCAA", "CCAC", "CCAG", "CCAT", "CCCA", "CCCC", "CCCG", "CCCT", 
        "CCGA", "CCGC", "CCGG", "CCGT", "CCTA", "CCTC", "CCTG", "CCTT", 
        "CGAA", "CGAC", "CGAG", "CGAT", "CGCA", "CGCC", "CGCG", "CGCT", 
        "CGGA", "CGGC", "CGGG", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT", 
        "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT", 
        "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT", 
        "GAAA", "GAAC", "GAAG", "GAAT", "GACA", "GACC", "GACG", "GACT", 
        "GAGA", "GAGC", "GAGG", "GAGT", "GATA", "GATC", "GATG", "GATT", 
        "GCAA", "GCAC", "GCAG", "GCAT", "GCCA", "GCCC", "GCCG", "GCCT", 
        "GCGA", "GCGC", "GCGG", "GCGT", "GCTA", "GCTC", "GCTG", "GCTT", 
        "GGAA", "GGAC", "GGAG", "GGAT", "GGCA", "GGCC", "GGCG", "GGCT", 
        "GGGA", "GGGC", "GGGG", "GGGT", "GGTA", "GGTC", "GGTG", "GGTT", 
        "GTAA", "GTAC", "GTAG", "GTAT", "GTCA", "GTCC", "GTCG", "GTCT", 
        "GTGA", "GTGC", "GTGG", "GTGT", "GTTA", "GTTC", "GTTG", "GTTT", 
        "TAAA", "TAAC", "TAAG", "TAAT", "TACA", "TACC", "TACG", "TACT", 
        "TAGA", "TAGC", "TAGG", "TAGT", "TATA", "TATC", "TATG", "TATT", 
        "TCAA", "TCAC", "TCAG", "TCAT", "TCCA", "TCCC", "TCCG", "TCCT", 
        "TCGA", "TCGC", "TCGG", "TCGT", "TCTA", "TCTC", "TCTG", "TCTT", 
        "TGAA", "TGAC", "TGAG", "TGAT", "TGCA", "TGCC", "TGCG", "TGCT", 
        "TGGA", "TGGC", "TGGG", "TGGT", "TGTA", "TGTC", "TGTG", "TGTT", 
        "TTAA", "TTAC", "TTAG", "TTAT", "TTCA", "TTCC", "TTCG", "TTCT", 
        "TTGA", "TTGC", "TTGG", "TTGT", "TTTA", "TTTC", "TTTG", "TTTT"]
    
    UTETRAMERS = dict.fromkeys(utetramers,0)
    for i in range(len(utetramers)):
        UTETRAMERS[utetramers[i]] = i
        
    return UTETRAMERS

def rcomp_kmers(x = get_utetramers().keys()):
    ''' Compute the reverse complement of kmers.
    
    @param  x a list of one or more kmers - assumed to be DNAAlphabet
    @return a dict of reverse_complement = kmer pairs
    '''
    
    from Bio.Seq import Seq
    from Bio.Alphabet import DNAAlphabet
    
    rc = dict()
    for k in x:
        rc[k] = str(Seq(k, DNAAlphabet()).reverse_complement())
    
    return rc