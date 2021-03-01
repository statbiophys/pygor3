#!/usr/bin/env python3
import pygor3 as p3


# def get_V_anchor_index_in_seq_from_gene(align:p3.IgorAlignment_data):
#     ins_size = len(align.insertions)
#     dels_size = len(align.deletions)
#     dict_gene_anchor_index = dict()
#     # TODO: COUNT number of insertions before cdr3_v_gene_anch


#
# # TODO: THERE IS A WAY TO GET
# def generate_str_fasta(indexed_sequence, list_vdj_alignments:dict):
#     """ Given an Sequence index and the corresponding alignments vj/ vdj
#     return a string with considering only offset"""
#
#     indexed_sequence.sequence = indexed_sequence.sequence.lower()
#     # add mismatches in sequence.
#     s = list(indexed_sequence.sequence)
#     for key_align in list_vdj_alignments.keys():
#         for pos_mis in list_vdj_alignments[key_align].mismatches:
#             s[pos_mis] = s[pos_mis].upper()
#     indexed_sequence.sequence = "".join(s)
#
#     str_fasta = ""
#     min_offset_key = min(list_vdj_alignments.keys(), key=lambda x: list_vdj_alignments[x].offset)  # .offset
#     min_offset = list_vdj_alignments[min_offset_key].offset
#     min_offset = min(indexed_sequence.offset, min_offset)
#
#     delta_offset = indexed_sequence.offset - min_offset
#     str_prefix = '-' * (delta_offset)
#     str_fasta_sequence = str_prefix + indexed_sequence.sequence
#     print(str_fasta_sequence)
#     str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
#     str_fasta = str_fasta + str_fasta_sequence + "\n"
#     for key in list_vdj_alignments.keys():
#         list_vdj_alignments[key].strGene_seq = list_vdj_alignments[key].strGene_seq.lower()
#         delta_offset = list_vdj_alignments[key].offset - min_offset
#         str_prefix = '-' * (delta_offset)
#         str_fasta_sequence = str_prefix + list_vdj_alignments[key].strGene_seq
#         print(str_fasta_sequence)
#         str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
#         str_fasta = str_fasta + str_fasta_sequence + "\n"
#         offset_5_p = list_vdj_alignments[key].offset_5_p - min_offset
#         offset_3_p = list_vdj_alignments[key].offset_3_p - min_offset
#         print("delta_offset : ", delta_offset)
#         print("offset_5_p : ", list_vdj_alignments[key].offset_5_p, offset_5_p)
#         print("offset_3_p : ", list_vdj_alignments[key].offset_3_p, offset_3_p)
#         str_prefix_2 = '-' * (offset_5_p+1)
#         str_fasta_sequence2 = str_prefix_2 + str_fasta_sequence[offset_5_p+1:offset_3_p+1]
#
#         str_fasta = str_fasta + "> " + list_vdj_alignments[key].strGene_name + ", score : "+str(list_vdj_alignments[key].score) + "\n"
#         str_fasta = str_fasta + str_fasta_sequence2 + "\n"
#
#         # TODO ADD MISMATCHES
#         align = list_vdj_alignments[key]
#         # align mismatches are in indexed sequence reference I need to convert it to gene reference given the alignment
#         # given the align.offset
#         # pos_in_gene  = pos_in_seq - align.offset
#         # pos_in_gene = cdr3 - align.offset
#
#     return str_fasta

def generate_csv_line(indexed_sequence:p3.IgorIndexedSequence, indexed_cdr3_record:list, list_vdj_alignments:dict, sep=';',
                      header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score']):
    csv_line = ""
    indexed_sequence.sequence = indexed_sequence.sequence.lower()
    fields_dict = dict()
    fields_dict['sequence_id'] = str(indexed_sequence.seq_index)
    fields_dict['sequence'] = indexed_sequence.sequence
    fields_dict['v_call'] =  list_vdj_alignments['V'].strGene_name
    fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
    fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
    fields_dict['v_score'] = str(list_vdj_alignments['V'].score)
    fields_dict['d_score'] = str(list_vdj_alignments['D'].score)
    fields_dict['j_score'] = str(list_vdj_alignments['J'].score)
    # FIXME: CREATE A BETTER WAY TO DO THIS
    fields_dict['v_anchor'] = ""
    fields_dict['j_anchor'] = ""
    fields_dict['junction'] = ""
    fields_dict['junction_aa'] = ""
    try:
        fields_dict['v_anchor'] = str(indexed_cdr3_record[1]) #""
        fields_dict['j_anchor'] = str(indexed_cdr3_record[2]) #""
        fields_dict['junction'] = str(indexed_cdr3_record[3]) #""
        fields_dict['junction_aa'] = str(indexed_cdr3_record[4]) #""
    except Exception as e:
        print("No junction for sequence : "+fields_dict['sequence_id'])
        print(e)
        pass

    # align = p3.IgorAlignment_data()

    for field in header_list:
        csv_line = csv_line + fields_dict[field]+sep
    return csv_line

    # mixcr :
    # targetSequences	targetQualities	allVHitsWithScore
    # allDHitsWithScore	allJHitsWithScore	allCHitsWithScore	allVAlignments	allDAlignments
    # allJAlignments	allCAlignments	nSeqFR1	minQualFR1	nSeqCDR1	minQualCDR1	nSeqFR2	minQualFR2	nSeqCDR2	minQualCDR2	nSeqFR3	minQualFR3	nSeqCDR3	minQualCDR3	nSeqFR4	minQualFR4	aaSeqFR1	aaSeqCDR1	aaSeqFR2	aaSeqCDR2	aaSeqFR3	aaSeqCDR3	aaSeqFR4	refPoints

    # AIRR :
    # sequence_id (req)
    # sequence (req)
    # rev_comp (req) //false by default # FIXME: for now
    # productive (req) // false by default # FIXME: a naive alignmnet of productive is possible.
    # locus (chain type) // eg. TCR IG
    # v_call
    # d_call
    # d2_call // FIXME: not gonna use this by the moment.
    # j_call
    # sequence alignment // FIXME: Aligned portion of query sequence. IMGT-gaps check the best way to do this
    # germline_alignment // FIXME: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any)
    # junction
    # junction_aa // Our CDR3
    # np1 (nts btwn V and D or btwn V and J  // vd_ins  or vj_ins
    # np2 (nts btwn first D and J gene or btwn first D and second D // dj_ins or d1d2_ins
    # np3 (nts btwn second D and J gene
    # v_score
    # v_cigar
    # d_score
    # d_cigar
    # j_score
    # j_cigar
    # junction_length




def generate_str_fasta_simple(indexed_sequence, list_vdj_alignments:dict):
    """ Given an Sequence index and the corresponding alignments vj/ vdj
    return a string with considering only offset"""

    str_fasta = ""
    str_fasta = str_fasta+"> "+ str(indexed_sequence.seq_index) + "\n"
    str_fasta = str_fasta+indexed_sequence.sequence+"\n"
    for key in list_vdj_alignments.keys():
        str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
        str_fasta = str_fasta + list_vdj_alignments[key].strGene_seq + "\n"

    # print(list_vdj_alignments['V'])
    # print(list_vdj_alignments['D'])
    return str_fasta



import argparse

def main():

    parser = argparse.ArgumentParser()
    # parser.add_argument("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    parser.add_argument("-o", "--output", dest="output", help='filename of csv file to export data')
    # parser.add_argument("-g", "--path_ref_genome", dest="path_ref_genome", help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv', default='./ref_genome')
    # parser.add_argument("-w", "--WD", dest="working_directory", help="Path where files gonna be created.", default='./')
    parser.add_argument("-D", "--database", dest="database", help="Igor database created with database script.")
    # parser.add_argument("")

    args = parser.parse_args()


    # # load gene templates
    # genomes = p3.IgorRefGenome()
    # genomes.path_ref_genome = args.path_ref_genome # "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/genomics_repseqio_F"
    # genomes.update_fln_names()
    # genomes.load_dataframes()


    # task = p3.IgorTask()
    # task.igor_wd = args.working_directory #"/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/test"
    # task.igor_batchname = args.batch #"sample"
    #
    # task.update_batch_filenames()

    db = p3.IgorSqliteDB()
    #db.flnIgorDB = task.igor_wd+"/"+task.igor_batchname+".db"
    db.fln_db = args.database
    db.connect_db()

    # Make a loop over all sequences
    if args.output is None:
        fln_output = args.database.split(".db")[0]+"_na.csv"
    else:
        fln_output = args.output

    ofile = open(fln_output, 'w')
    seq_index_list = db.execute_select_query("SELECT seq_index FROM IgorIndexedSeq;")
    seq_index_list = map(lambda x:x[0], seq_index_list)
    header_list = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score', 'v_anchor', 'j_anchor', 'junction', 'junction_aa']
    sep = ";"
    str_header = sep.join(header_list)
    ofile.write(str_header+'\n')
    for seq_index in seq_index_list:
        try:
            # seq_index = args.seq_index
            indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
            indexed_sequence.offset = 0

            indexed_cdr3_record = db.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
            print(indexed_cdr3_record)

            best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
            best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
            best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)

            vdj_naive_alignment = {'V': best_v_align_data,
                               'D': best_d_align_data,
                               'J': best_j_align_data}

            v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
            # print('V', len(v_align_data_list), [ ii.score for ii in v_align_data_list])
            try:
                d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
            except Exception as e:
                print(e)
                print("No D sequences alignments found!")
                pass
            # print('D', len(d_align_data_list), [ ii.score for ii in d_align_data_list])
            j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
            # print('J', len(j_align_data_list), [ ii.score for ii in j_align_data_list])

            # 1. Choose the highest score then check if this one is the desire range.
            # if there is an overlap
            # calculate score without overlap. If overlap
            # if hightest score
            for i, d_align_data in enumerate(d_align_data_list):
                # Check if D is btwn V and J position
                if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (d_align_data.offset_3_p <= best_j_align_data.offset_5_p) :
                    #vdj_naive_alignment['D'+str(i)] = d_align_data
                    vdj_naive_alignment['D'] = d_align_data
                    break

            # ofile = open(task.igor_batchname+'__'+str(indexed_sequence.seq_index)+'_na.fasta', 'w')
            # str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
            # ofile.write(str_fasta)

            str_csv_line = generate_csv_line(indexed_sequence, indexed_cdr3_record, vdj_naive_alignment, sep=sep, header_list=header_list)
            ofile.write(str_csv_line+"\n")
        except Exception as e:
            print("balsa : "+str(seq_index))
            print(e)
            pass
    ofile.close()

    # given a IgorGenomic template and using pygor database for IGoR get anchors from alignment and then

    # algn_fln = 'sample__2_na.fasta'
    # seq_lens = [len(record.seq) for record in SeqIO.parse(algn_fln, 'fasta')]
    # max_seq_lens = max(seq_lens)
    # new_alignments = list()
    # for record in SeqIO.parse(algn_fln, 'fasta'):
    #     str_extra = max_seq_lens - len(record.seq)
    #     record.seq = record.seq + "-" * str_extra
    #     new_alignments.append(record)
    #     # print(len(record.seq), max_seq_lens)
    #
    # tmp_fln = 'example.fasta'
    # SeqIO.write(new_alignments, tmp_fln, "fasta")

    # da_heavy_pen_nuc44 = p3.utils.da_heavy_pen_nuc44_vect
    # print(da_heavy_pen_nuc44)
    # import pandas as pd
    # df = pd.DataFrame(data=da_heavy_pen_nuc44.values, index=da_heavy_pen_nuc44['lbl__' + 'x'].values,
    #                   columns=da_heavy_pen_nuc44['lbl__' + 'y'].values)
    #
    # print(df)


if __name__ == "__main__":
    main()
