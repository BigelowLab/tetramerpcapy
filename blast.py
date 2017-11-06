import subprocess
import sys
from distutils.spawn import find_executable
from xml.etree import cElementTree as ElementTree


def which_blast(app = 'blastn'):
    '''
    Find an application by name
    
    @param      app the name of the application to test
    @return     the fully qualified path or None
    '''
    return find_executable(app)
    
    
    
def has_blast(app = 'blastn'):
    '''
    Determine is an application (executable) is available
    
    @param      app the name of the application to test
    @return     boolean, True if the application is found
    '''
    return which_blast(app) is not None
    

def run_blast(cmd = 'blastn -db nt -num_threads 7 -num_alignments 10 -evalue 10',
    in_file = 'example-outliers.fasta',
    out_file = 'example-outliers.xml',
    verbose = False):
    '''
    Run the blast executable
    
    @param cmd          the command srting less input/output info
    @param in_file      the name of the input fasta file
    @param out_file     the name of the xml_file to output
    @param verbose      boolean, if True then echo the blast command string
    @return             the value returned by subprocess.call           
    '''
             
    arglist = cmd.split(' ') + ['-outfmt 5', '-query',
        in_file, '-out',  out_file]
    
    argstring = ' '.join(arglist)
    
    if verbose:
        print('CMD: ', argstring)
        

    ok = subprocess.call(argstring, shell = True)
    return ok
    

def blast_to_table(in_file = 'example-outliers.xml',
    out_file = 'example-outliers.tsv', sep = '\t'):
    
    '''
    Dump a BLAST XML file to tabular output.
    
    Based upon blastxml_to_tabular.py. No processing occurs in the dump process.
    Column names, found in the first row, are modified from the tag-names in the 
    Blast XML output nodes: dashes, '-', are replaced with underscores, '_'.
      Iteration_iter_num
      Iteration_query_id
      Iteration_query_def
      Iteration_query_len
      Hit_num
      Hit_id
      Hit_def
      Hit_accession
      Hit_len
      Hsp_num
      Hsp_bit_score
      Hsp_score
      Hsp_evalue
      Hsp_query_from
      Hsp_query_to
      Hsp_hit_from
      Hsp_hit_to
      Hsp_query_frame
      Hsp_identity
      Hsp_positive
      Hsp_gaps
      Hsp_align_len
      Hsp_qseq
      Hsp_hseq
      Hsp_midline
      
    @param  in_file   input XML filename
    @param  out_file  output CSV filename
    @param  integer number of records saved
    '''
    
    def stop_err( msg ):
      sys.stderr.write('%s\n' % msg)
      return 0
    
    
    # get an iterable
    try: 
      context = ElementTree.iterparse(in_file, events=('start', 'end'))
    except:
      stop_err('Invalid input file data format.')
    # turn it into an iterator
    context = iter(context)
    # get the root element
    try:
      event, root = context.next()
    except:
      stop_err( 'Invalid input file data format.' )
    
    # write a nice header to the output file, leave the ouptut file open for now
    header = ['Iteration_iter_num',
              'Iteration_query_id',
              'Iteration_query_def',
              'Iteration_query_len',
              'Hit_num',
              'Hit_id',
              'Hit_def',
              'Hit_accession',
              'Hit_len',
              'Hsp_num',
              'Hsp_bit_score',
              'Hsp_score',
              'Hsp_evalue',
              'Hsp_query_from',
              'Hsp_query_to',
              'Hsp_hit_from',
              'Hsp_hit_to',
              'Hsp_query_frame',
              'Hsp_identity',
              'Hsp_positive',
              'Hsp_gaps',
              'Hsp_align_len',
              'Hsp_qseq',
              'Hsp_hseq',
              'Hsp_midline',]
    outfile = open(out_file, 'w')
    outfile.write(sep.join(header) + '\n')
    blast_program = None
    for event, elem in context:
      if event == 'end' and elem.tag == 'BlastOutput_program':
        blast_program = elem.text
      # for every <Iteration> tag
      if event == 'end' and elem.tag == 'Iteration':
        #Expecting either this, from BLAST 2.2.25+ using FASTA vs FASTA
        # <Iteration_query-ID>sp|Q9BS26|ERP44_HUMAN</Iteration_query-ID>
        # <Iteration_query-def>Endoplasmic reticulum resident protein 44 OS=Homo sapiens GN=ERP44 PE=1 SV=1</Iteration_query-def>
        # <Iteration_query-len>406</Iteration_query-len>
        # <Iteration_hits></Iteration_hits>
        #
        #Or, from BLAST 2.2.24+ run online
        # <Iteration_query-ID>Query_1</Iteration_query-ID>
        # <Iteration_query-def>Sample</Iteration_query-def>
        # <Iteration_query-len>516</Iteration_query-len>
        # <Iteration_hits>...
        
        iter_num = elem.findtext('Iteration_iter-num')
        query_ID = elem.findtext('Iteration_query-ID')
        query_def = elem.findtext('Iteration_query-def')
        query_len = elem.findtext('Iteration_query-len')
    
        # for every <Hit> within <Iteration>
        for hit in elem.findall('Iteration_hits/Hit'):
          #Expecting either this,
          # <Hit_id>gi|3024260|sp|P56514.1|OPSD_BUFBU</Hit_id>
          # <Hit_def>RecName: Full=Rhodopsin</Hit_def>
          # <Hit_accession>P56514</Hit_accession>
          #or,
          # <Hit_id>Subject_1</Hit_id>
          # <Hit_def>gi|57163783|ref|NP_001009242.1| rhodopsin [Felis catus]</Hit_def>
          # <Hit_accession>Subject_1</Hit_accession>
          #
          #apparently depending on the parse_deflines switch
          #
          #Or, with BLAST 2.2.28+ can get this,
          # <Hit_id>gnl|BL_ORD_ID|2</Hit_id>
          # <Hit_def>chrIII gi|240255695|ref|NC_003074.8| Arabidopsis thaliana chromosome 3, complete sequence</Hit_def>
          # <Hit_accession>2</Hit_accession>
          hit_num = hit.findtext('Hit_num')
          hit_ID = hit.findtext('Hit_id')
          hit_def = hit.findtext('Hit_def')
          hit_accession = hit.findtext('Hit_accession')
          hit_len = hit.findtext('Hit_len')
          
          # for every <Hsp> within <Hit>
          for hsp in hit.findall('Hit_hsps/Hsp'):
            values = [iter_num, 
                      query_ID, 
                      query_def, 
                      query_len,
                      hit_num, 
                      hit_ID, 
                      hit_def, 
                      hit_accession, 
                      hit_len,
                      hsp.findtext('Hsp_num'),
                      hsp.findtext('Hsp_bit-score'),
                      hsp.findtext('Hsp_score'),
                      hsp.findtext('Hsp_evalue'),
                      hsp.findtext('Hsp_query-from'),
                      hsp.findtext('Hsp_query-to'),
                      hsp.findtext('Hsp_hit-from'),
                      hsp.findtext('Hsp_hit-to'),
                      hsp.findtext('Hsp_query-frame'),
                      hsp.findtext('Hsp_identity'),
                      hsp.findtext('Hsp_positive'),
                      hsp.findtext('Hsp_gaps'),
                      hsp.findtext('Hsp_align-len'),
                      hsp.findtext('Hsp_qseq'),
                      hsp.findtext('Hsp_hseq'),
                      hsp.findtext('Hsp_midline'),] 
            outfile.write(sep.join(values) + '\n')
        # prevents ElementTree from growing large datastructure
        root.clear()
        elem.clear()
    outfile.close()
    
    return len(values)
    
    