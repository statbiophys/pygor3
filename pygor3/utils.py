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