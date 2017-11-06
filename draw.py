import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import re
import misc
import datetime
from adjustText import adjust_text


def get_explanation(X):
    '''
    Retrieve the normalized 'explanation' of variance for each PC as eigenvalues
    
    Returns the 'explanations' of variance as derived from the fracs component of the 
     PC (principal components)
    
    @param X the tetramer dict to rule them all
    @return The proportion of variance of each of the principal components formatted as string.
    '''
    p = list(X['P'].fracs)
    
    r = dict()
    for i in range(1, X['params']['npc']+1):
        r['PC{0}'.format(i)] = ', %0.1f%% variation explained' % (p[i] * 100)
    
    return r
    
def get_tetramer_range(X):
    ''' 
    Compute the min,max ranges for each PC loading

    @param  X the tetra dist
    @return pandas dataframe with one row per PC and two columns (min and max)
    '''
    mx = X['x'].max(0)
    mx.name = 'max'
    mn = X['x'].min(0)
    mn.name = 'min'
    
    return pd.concat([mn, mx], axis = 1)

def trim_legend_text(x = 'foo\t bar baz tom petty had\n a pretty good run',
    max_len = 35, strip_white = True, ext = ' ...'):
    '''
    Convert the provided string to a shortened prettified verision for graphics

    @param  max_len, the maximum output string length
    @param  strip_white, removes tabs and newlines
    @param  ext, append this when input string has to be trimmed
    @return the input string possibly trimmed and prettified
    '''
    
    pat = re.compile(r'[\t\n]')
    x = re.sub(pat, '', x)
    
    if (len(x) >= max_len):
        x = x[1:(max_len - len(ext))] + ext
    
    return x

def blast_hit_text(x, max_len = 30, no_hit_text = "no significant hits"):
    '''
    Construct strings suitable for drawing
    
    hit_accession blast score = n
    hit_def_up_to_30_characters
    
    @param x a dataframe of the blast data
    @param max_len the maximum allowed characters for the second row of text
    @param no_hit_text substituted text fro when hit_len <= 0
    @return dict of windowname: text_to_draw, one element per row of input dataframe
    '''
    
    stub = "%0.10s blast score = %0.0f\n%0." + str(max_len) + "s"
    def do_one(x, stub = '%0.10s blast score = %0.0f\n%0.30s'):
        txt = stub % (x['Hit_accession'], float(x['Hsp_bit_score']),x['Hit_def'] )
        if x['Hit_len'] <= 0:
            txt = no_hit_text
        return txt
    
    txt = [do_one(x.iloc[i]) for i in range(0, x.shape[0])]
    nm = list(x.index.values)
    dct = dict(zip(nm, txt))
    return dct
    
    
def plot_PCvPC(X, iplot = 0, pdf = None, add_blast = True):
    
    colors = {
        'gray': "#BFBFBF", 'blue': '#0000FF',
        'py': ["#D01C8B", "#F1B6DA", "#B8E186", "#4DAC26"],
        'px': ["#E66101", "#FDB863", "#B2ABD2", "#5E3C99"]
        }
    markers = {
        'gray': 'o',
        'small': 4,
        'big':   6,
        'px':  ["^", "^", "v", "v"], 
        'py':  ["^", "^", "v", "v"]
        }        
    lim = X['plot']['lim']
    
    # make an array ala [PC1, PC2], [PC3, PC4], ...
    nmpc = X['plot']['nmpc']
    
    nplots = nmpc.shape[0]
    blastinfo = X['params']['blastcmd'].split(" ")
    date = datetime.date.today().strftime("%Y-%m-%d")
    subtitle_stub = 'window = {0}  step = {1}\n{2} against GenBank {3}, {4}'
    
    out = X['plot']['out'] 
    x = X['plot']['x']    
    
    fig = plt.figure(figsize=(8, 8))
    
    # these yield 'PC1', 'PC2'  etc
    px = nmpc[iplot,0]
    py = nmpc[iplot,1]
    
    # define a subplot so we can shrink it later
    ax = plt.subplot(111,aspect = 'auto')
        
    # these are just the axis labels
    expl = get_explanation(X)
    
    plt.plot(x[px], x[py], 'o',
        marker = markers['gray'], markersize = markers['small'], 
        markerfacecolor = colors['gray'], 
        markeredgecolor = colors['gray'])
    
            
    plt.suptitle(X['name'], y= 0.98, fontsize=18)
    plt.title(subtitle_stub.format(X['params']['window'], X['params']['step'], blastinfo[0], blastinfo[2], date),
        fontsize=10)
    plt.xlabel(''.join([px,expl[px]]))
    plt.ylabel(''.join([py,expl[py]]))
    
    # shrink the plotting box to make room for the legend below
    dy = 0.8
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * (1-dy),
                     box.width, box.height * dy])
    
    xpad = [-10, 10]
    xlim = ax.get_xlim()
    new_xlim = [xlim[0] + xpad[0], xlim[1] + xpad[1]]
    ax.set_xlim(new_xlim)
    
    print("old" , '{} {}'.format(xlim[0], xlim[1]))
    print("new" , '{} {}'.format(new_xlim[0], new_xlim[1]))
    
    
    # now get ready to plot outliers and capture legend info
    # the references to x and y are a bit flummoxed here.  We select outliers in each dimension (each PC)
    # and it is these that are refered to x and y. In hindsight, I wish I called them One and Two (or A and B.)
    # Anyway, each outliers can be plotted in x and y regardless if it is One or Two
    
    # dataframes for One and Two
    outx = out[out.PC == px]
    outy = out[out.PC == py]
    
    #contig names for One and Two
    contigx = list(outx.index.values)
    contigy = list(outy.index.values)
    
    # fudge the location of the column titles for One and Two
    dystep = 0.03
    toptext = 1 - dy - dystep
    lefttext = .12  # these are guesses
    righttext = .63 # for providing a dual title to the legend
    
    plt.annotate(px, 
                 xy=(lefttext, toptext),
                 xytext=(lefttext, toptext),  
                 xycoords='figure fraction',
                 textcoords='figure fraction',
                 arrowprops=None)
    plt.annotate(py, 
                 xy=(righttext, toptext),
                 xytext=(righttext, toptext),  
                 xycoords='figure fraction',
                 textcoords='figure fraction',
                 arrowprops=None)  
      
    # a list for legend elements [One.1, One.2. One.3. One.4, Two.1, Two.2, Two.3, Two.4]
    leg = [] 
    
    # plot points for One         
    for j in range(0, len(contigx)):
        
        a, = ax.plot(x.loc[contigx[j]][px], x.loc[contigx[j]][py],
            marker = markers['px'][j], markersize = markers['big'],
            color = colors['px'][j],
            markeredgecolor = colors['px'][j], 
            markerfacecolor = colors['px'][j],
            linestyle = '-',
            label = contigx[j]) 
        leg.append(a)
     
    # plot points for Two   
    for j in range(0, len(contigy)):
        a, = ax.plot(x.loc[contigy[j]][px], x.loc[contigy[j]][py],
            marker = markers['py'][j], markersize = markers['big'],
            color = colors['py'][j],
            markeredgecolor = colors['py'][j],
            markerfacecolor = colors['py'][j],
            linestyle = '-')
        leg.append(a)        
    
    # make the legend and colorize the text
    leg = ax.legend(leg, contigx + contigy,
        loc='lower center', bbox_to_anchor = (0.1, 0., 0.8, .05),
        bbox_transform = plt.gcf().transFigure,
        mode="expand", borderaxespad=0.,
        fancybox=False, shadow=False, ncol=2,
        numpoints = 1, frameon = False)
    for color,text in zip(colors['px'] + colors['py'],leg.get_texts()):
        text.set_color(color)
      
      
    if add_blast: 
        # get a dict of the blast text 
        blast_text = blast_hit_text(X['xblast'])
        y = x.copy()
        y.index = y.index.droplevel('name')
        texts = dict()
        for k in blast_text.keys():
            ha = 'left' if y.loc[k,px] > 0 else 'right'
            txt = ax.text(y.loc[k,px], y.loc[k,py], blast_text[k],
                color = colors['blue'], size = 'x-small',
                horizontalalignment = ha)
            texts[k] = txt
        adjust_text(texts.values())
        
    pdf.savefig()
    plt.close()       
    return None
    
def plot_tetramer(X, filename = 'example-tetramer.pdf', add_blast = True):
    '''
    Plot a 4-page graphic that shows PC1 v PC2, PC3 v PC4, PC5 v PC6, PC7 v PC8 
    with outliers highlghted and labeled with bast results
    
    @seealso http://matplotlib.org/1.5.1/faq/howto_faq.html#save-multiple-plots-to-one-pdf-file
    @seealso http://matplotlib.org/examples/pylab_examples/multipage_pdf.html
    
    @param X        dict of tetramer analysis results
    @param out_file name of the output PDF file
    '''
 
            
 
    lim = get_tetramer_range(X)
    # make an array ala [PC1, PC2], [PC3, PC4], ...
    nmpc = lim.index.values.reshape(X['params']['npc']/2,2) 
    nplots = nmpc.shape[0]
    
    out = X['outliers'].copy()
    
    x = X['x'].copy()
    ix = misc.parse_windowname(list(x.index.values))
    ix1 = list(ix['name'].values)
    ix2 = list(ix.index.values)
    idx =  zip(ix1,ix2)
    x['name'] = ix1
    x['window'] = ix2
    x = x.set_index(['name','window'])
    
    
    X['plot'] = {
        'lim'   : lim,
        'nmpc'  : nmpc,
        'nplots': nplots,
        'out'   : out,
        'x'     : x
    }
    
    
    
    with PdfPages(filename) as pdf:
        
        for iplot in range(0, 4):
            plot_PCvPC(X, iplot = iplot, pdf = pdf, add_blast = add_blast)
