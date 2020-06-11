#!/usr/bin/env python3
import pygor3 as p3


def mismatches_upper(str_sequence, mismatches):
    s = list(str_sequence.lower())
    for pos_mis in mismatches:
        s[pos_mis] = s[pos_mis].upper()
    str_sequence_mis_upper = "".join(s)
    return str_sequence_mis_upper

from  optparse import OptionParser
def main():
    #parser = OptionParser()
    #parser.add_option("-s", "--species", dest="species", help='Igor species')
    #parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    #parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    #parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    #(options, args) = parser.parse_args()
    #print("hola fasdasdfd 888888")
    # species="human"
    # chain="tcr_beta"
    # mdl = p3.IgorModel.load_default(species, chain)
    # import hvplot
    # hvplot.save(mdl.parms.plot_Graph(), species+"__"+chain+".html")

    mdl = p3.IgorModel.load_default("human", "tcr_beta")
    print( mdl.get_events_nicknames_list() )
    import networkx as nx
    print( list( nx.topological_sort(mdl.parms.G) ))
    print("SORTED")
    print( list( map(lambda x: x.nickname, mdl.parms.get_Event_list_sorted()) ) )
    # xdata_list_roots = list()
    # for event_nickname in mdl.xdata.keys():
    #     if len(mdl.xdata[event_nickname].parents) == 0 :
    #         xdata_list_roots.append(mdl.xdata[event_nickname])
    #         #xdata_list.append(mdl.parms.get_Event(event_nickname))
    # print(xdata_list_roots)
    # print(list(map(lambda x: x.nickname, roots)))
    # roots_sorted = sorted(roots, key=lambda x: x.nickname)
    # print(roots_sorted)
    # print(list(map(lambda x: x.nickname, roots_sorted)))
    # ordered_list = list()
    # for root in roots:
    #     # first process the parents before the childs.
    #     ordered_list.append(root)

    mdl.parms.Event_list



"""
    # load gene templates
    genomes = p3.IgorRefGenome()
    genomes.path_ref_genome = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/genomics_repseqio_F"
    genomes.update_fln_names()
    genomes.load_dataframes()

    # print(genomes.df_V_ref_genome['anchor_index'])

    db = p3.IgorSqliteDB()
    # db.flnIgorDB = task.igor_wd + "/" + task.igor_batchname + ".db"
    igor_wd = "test"
    igor_batchname = "sample"
    db.flnIgorDB = igor_wd + "/" + igor_batchname + ".db"
    db.connect_db()

    indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(0)# p3.IgorIndexedSequence.
    v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
    # print('-'*50)

    v_align_min_offset = min(v_align_data_list, key=lambda x: x.offset)  # .offset
    print('='*50)
    # print(v_align_min_offset)
    min_offset = v_align_min_offset.offset

    df = genomes.df_V_ref_genome.set_index('name')
    dict_anchor_index = genomes.get_anchors_dict()
    for v_align in v_align_data_list:
        print('*' * 20)
        # 1. get anchor position in sequence from gene template
        anchor_index_in_gene = dict_anchor_index['V'][v_align.strGene_name]
        anchor_index_in_seq = anchor_index_in_gene + v_align.offset  # + insertions or deletions before. that position.

        indexed_sequence.sequence = indexed_sequence.sequence.lower()

        # Show
        s = list(indexed_sequence.sequence.lower())
        for pos_mis in v_align.mismatches:
            s[pos_mis] = s[pos_mis].upper()
        str_sequence_fasta = "".join(s)

        delta_offset = 0 - min_offset
        str_prefix = '-' * (delta_offset)
        str_sequence_fasta = str_prefix + str_sequence_fasta
        print(str_sequence_fasta)
        delta_offset = v_align.offset - min_offset
        str_prefix = ' ' * (delta_offset)
        str_sequence_fasta = str_prefix + v_align.strGene_seq.lower()
        print(str_sequence_fasta)

    print(v_align)
    # task = p3.IgorTask()
    # task.igor_wd = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/test"
    # task.igor_batchname = "sample"
    #
    # task.update_batch_filenames()
    #
    # db = p3.IgorSqliteDB()
    
"""

    # load Anchors with anchors extract CDR3 and naive alignment

if __name__ == "__main__":
    main()
