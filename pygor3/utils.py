import collections
def from_fasta_to_dataframe(fln_fasta):
    # Fasta to dataframe
    from Bio import SeqIO
    import pandas as pd
    genes_name_list = list()
    genes_value_list = list()
    for gene_record in SeqIO.parse(fln_fasta, "fasta"):
        genes_name_list.append(gene_record.description)
        genes_value_list.append(str(gene_record.seq))
    df_genes = pd.DataFrame.from_dict({'name': genes_name_list, 'value': genes_value_list})
    df_genes.index.name = 'id'
    return df_genes



# // A, C, G, T, R, Y, K, M, S, W, B, D, H, V, N
heavy_pen_nuc44_vect = [
5, -14, -14, -14, -14, 2, -14, 2, 2, -14, -14, 1, 1, 1, 0,
-14, 5, -14, -14, -14, 2, 2, -14, -14, 2, 1, -14, 1, 1, 0,
-14, -14, 5, -14, 2, -14, 2, -14, 2, -14, 1, 1, -14, 1, 0,
-14, -14, -14, 5, 2, -14, -14, 2, -14, 2, 1, 1, 1, -14, 0,
-14, -14, 2, 2, 1.5, -14, -12, -12, -12, -12, 1, 1, -13, -13, 0,
2, 2, -14, -14, -14, 1.5, -12, -12, -12, -12, -13, -13, 1, 1, 0,
-14, 2, 2, -14, -12, -12, 1.5, -14, -12, -12, 1, -13, -13, 1, 0,
2, -14, -14, 2, -12, -12, -14, 1.5, -12, -12, -13, 1, 1, -13, 0,
2, -14, 2, -14, -12, -12, -12, -12, 1.5, -14, -13, 1, -13, 1, 0,
-14, 2, -14, 2, -12, -12, -12, -12, -14, 1.5, 1, -13, 1, -13, 0,
-14, 1, 1, 1, 1, -13, 1, -13, -13, 1, 0.5, -12, -12, -12, 0,
1, -14, 1, 1, 1, -13, -13, 1, 1, -13, -12, 0.5, -12, -12, 0,
1, 1, -14, 1, -13, 1, -13, 1, -13, 1, -12, -12, 0.5, -12, 0,
1, 1, 1, -14, -13, 1, 1, -13, 1, -13, -12, -12, -12, 0.5, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

import xarray as xr
import numpy as np
list_nt_lbl = ['A', 'C', 'G', 'T', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']
da_heavy_pen_nuc44_vect = xr.DataArray(np.array(heavy_pen_nuc44_vect).reshape(15, 15), \
                               dims=('x', 'y'))
# print(len(list_nt_lbl))
strDim = 'x'
da_heavy_pen_nuc44_vect[strDim] = range(len(list_nt_lbl))
strCoord = 'lbl__' + strDim
da_heavy_pen_nuc44_vect[strCoord] = (strDim, list_nt_lbl)

strDim = 'y'
da_heavy_pen_nuc44_vect[strDim] = range(len(list_nt_lbl))
strCoord = 'lbl__' + strDim
da_heavy_pen_nuc44_vect[strCoord] = (strDim, list_nt_lbl)



#####################################################
########## PLOTTING SEQUENCES  #################
#####################################################
# Functions copied from :
# Reference: https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner

# FIXME: PLEASE MAKE IT MORE MATPLOTLIB-ISH
try:
    from Bio import AlignIO, SeqIO
    import numpy as np
    from bokeh.plotting import figure, output_file, save, show
    from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
    from bokeh.models.glyphs import Text, Rect
    from bokeh.layouts import gridplot

    def get_colors(seqs):
        """make colors for bases in sequence"""
        text = [i for s in list(seqs) for i in s]
        clrs = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', '-': 'white',
                'a': 'red', 't': 'green', 'g': 'orange', 'c': 'blue'}
        colors = [clrs[i] for i in text]
        return colors


    def view_alignment(aln, fontsize="9pt", plot_width=800):
        """Bokeh sequence alignment view"""
        # make sequence and id lists from the aln object
        seqs = [rec.seq for rec in (aln)]
        ids = [rec.id for rec in aln]
        text = [i for s in list(seqs) for i in s]

        colors = get_colors(seqs)
        N = len(seqs[0])
        S = len(seqs)
        width = .4
        x = np.arange(1, N + 1)
        y = np.arange(0, S, 1)
        # creates a 2D grid of coords from the 1D arrays
        xx, yy = np.meshgrid(x, y)
        # flattens the arrays
        gx = xx.ravel()
        gy = yy.flatten()
        # use recty for rect coords with an offset
        recty = gy + .5
        h = 1 / S
        # now we can create the ColumnDataSource with all the arrays
        source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
        plot_height = len(seqs) * 15 + 50
        x_range = Range1d(0, N + 1, bounds='auto')
        if N > 100:
            viewlen = 100
        else:
            viewlen = N
        # view_range is for the close up view
        view_range = (0, viewlen)
        tools = "xpan, xwheel_zoom, reset, save"

        # entire sequence view (no text, with zoom)
        p = figure(title=None, plot_width=plot_width, plot_height=50,
                   x_range=x_range, y_range=(0, S), tools=tools,
                   min_border=0, toolbar_location='below')
        rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                     line_color=None, fill_alpha=0.6)
        p.add_glyph(source, rects)
        p.yaxis.visible = False
        p.grid.visible = False

        # sequence text view with ability to scroll along x axis
        p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
                    x_range=view_range, y_range=ids, tools="xpan,reset",
                    min_border=0, toolbar_location='below')  # , lod_factor=1)
        glyph = Text(x="x", y="y", text="text", text_align='center', text_color="black",
                     text_font="monospace", text_font_size=fontsize)
        rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                     line_color=None, fill_alpha=0.4)
        p1.add_glyph(source, glyph)
        p1.add_glyph(source, rects)

        p1.grid.visible = False
        p1.xaxis.major_label_text_font_style = "bold"
        p1.yaxis.minor_tick_line_width = 0
        p1.yaxis.major_tick_line_width = 0

        p = gridplot([[p], [p1]], toolbar_location='below')
        show(p)
        return p


except ImportError as error:
    # Output expected ImportErrors.
    # print(error.__class__.__name__ + ": " + error.message)
    # FIXME: no error printing for autocomplete
    # print(error)
    pass
except Exception as exception:
    # Output unexpected Exceptions.
    print(exception, False)
    print(exception.__class__.__name__ + ": " + exception.message)


def run_igor_datadir():
    import subprocess

    cmd = "which igor"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    igor_exec_path = line.decode("utf-8").replace('\n', '')

    cmd = igor_exec_path + " -getdatadir"

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = []
    while True:
        line = p.stdout.readline()
        line = line.decode("utf-8")
        stdout.append(line)
        # print (line, end='')
        if line == '' and p.poll() != None:
            break
    return (''.join(stdout)).replace('\n','')


class GeneSegment:
    def __int__(self, gene_type=None):
        self.gene_type = gene_type
        self.palindrome_5_end = None
        self.gene_ini = None
        self.gene_end = None
        self.gene_cut = None
        self.palindrome_3_end = None
        self.gene_segment = None

class InsertSegment:
    def __int__(self, gene_type=None):
        self.gene_type = gene_type
        self.palindrome_5_end = None
        self.gene_ini = None
        self.gene_end = None
        self.gene_cut = None
        self.palindrome_3_end = None
        self.gene_segment = None

def get_gene_segment(str_gene_template, int_gene_5_del=None, int_gene_3_del=None):
    if int_gene_5_del is None:
        int_gene_5_del = 0
    if int_gene_3_del is None:
        int_gene_3_del = 0

    # print(str_gene_template, int_gene_5_del, int_gene_3_del)

    int_ini = 0
    int_end = len(str_gene_template)
    str_gene_3_palindrome = ""
    str_gene_5_palindrome = ""

    if int_gene_5_del < 0:
        int_ini = 0
        str_gene_5_palindrome = dna_complementary( (str_gene_template[:-int_gene_5_del])[::-1] )
    else:
        int_ini = int_gene_5_del
        str_gene_5_palindrome = ""

    if int_gene_3_del < 0:
        int_end = len(str_gene_template)
        str_gene_3_palindrome = dna_complementary( (str_gene_template[int_gene_3_del:])[::-1] )
    else:
        int_end = len(str_gene_template) - int_gene_3_del
        str_gene_3_palindrome = ""

    segment_dict = collections.OrderedDict()
    segment_dict['palindrome_5_end'] = str_gene_5_palindrome
    segment_dict['gene_ini'] = int_ini
    segment_dict['gene_end'] = int_end
    segment_dict['gene_cut'] = str_gene_template[int_ini:int_end]
    segment_dict['palindrome_3_end'] = str_gene_3_palindrome
    segment_dict['gene_segment'] = str_gene_5_palindrome + str_gene_template[int_ini:int_end] + str_gene_3_palindrome
    return segment_dict
    # str_gene_segment = str_gene_5_palindrome + str_gene_template[int_ini:int_end] + str_gene_3_palindrome
    # return str_gene_segment


def dna_complementary(str_seq):
    from Bio.Seq import Seq
    return str(Seq(str_seq).complement())
