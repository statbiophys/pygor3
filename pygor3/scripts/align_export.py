#!/usr/bin/env python3
import pygor3 as p3


# TODO: THERE IS A WAY TO GET
def generate_str_fasta(indexed_sequence, list_vdj_alignments:dict):
    """ Given an Sequence index and the corresponding alignments vj/ vdj
    return a string with considering only offset"""

    str_fasta = ""
    min_offset_key = min(list_vdj_alignments.keys(), key=lambda x: list_vdj_alignments[x].offset)  # .offset
    min_offset = list_vdj_alignments[min_offset_key].offset
    min_offset = min(indexed_sequence.offset, min_offset)

    delta_offset = indexed_sequence.offset - min_offset
    str_prefix = '-' * (delta_offset)
    str_fasta_sequence = str_prefix + indexed_sequence.sequence
    print(str_fasta_sequence)
    str_fasta = str_fasta+"> "+ str(indexed_sequence.seq_index) + "\n"
    str_fasta = str_fasta+str_fasta_sequence+"\n"
    for key in list_vdj_alignments.keys():
        delta_offset = list_vdj_alignments[key].offset - min_offset
        str_prefix = '-' * (delta_offset)
        str_fasta_sequence = str_prefix + list_vdj_alignments[key].strGene_seq
        print(str_fasta_sequence)
        str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
        str_fasta = str_fasta + str_fasta_sequence + "\n"
        offset_5_p = list_vdj_alignments[key].offset_5_p - min_offset
        offset_3_p = list_vdj_alignments[key].offset_3_p - min_offset
        print("delta_offset : ", delta_offset)
        print("offset_5_p : ", list_vdj_alignments[key].offset_5_p, offset_5_p)
        print("offset_3_p : ", list_vdj_alignments[key].offset_3_p, offset_3_p)
        str_prefix_2 = '-' * (offset_5_p)
        str_fasta_sequence2 = str_prefix_2 + str_fasta_sequence[offset_5_p:offset_3_p]

        str_fasta = str_fasta + "> " + list_vdj_alignments[key].strGene_name + "\n"
        str_fasta = str_fasta + str_fasta_sequence2 + "\n"


        # aaa = seq_tmp[(item.offset_5_p - min_offset):(item.offset_3_p - min_offset) + 1]
        # print(item.strGene_class, item.offset_5_p, item.offset_3_p)
        # print(aaa)

    # print(list_vdj_alignments['V'])
    # print(list_vdj_alignments['D'])
    return str_fasta


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


from  optparse import OptionParser
def main():
    parser = OptionParser()
    parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    # load gene templates
    genomes = p3.IgorRefGenome()
    genomes.path_ref_genome = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/genomics_repseqio_F"
    genomes.update_fln_names()
    df_V, df_J = genomes.load_dataframes()

    task = p3.IgorTask()
    task.igor_wd = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/test"
    task.igor_batchname = "sample"

    task.update_batch_filenames()

    db = p3.IgorSqliteDB()
    flnIgorDB = task.igor_batchname+".db"
    db.flnIgorSQL = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/pygor3/IgorDB.sql"
    db.createSqliteDB(flnIgorDB)
    # db.load_VDJ_Database(
    #     task.igor_fln_indexed_sequences,
    #     genomes.fln_genomicVs,
    #     genomes.fln_genomicDs,
    #     genomes.fln_genomicJs,
    #     task.igor_fln_align_V_alignments,
    #     task.igor_fln_align_D_alignments,
    #     task.igor_fln_align_J_alignments)

    seq_index = 0
    indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
    indexed_sequence.offset = 0

    best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
    best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
    best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)

    vdj_naive_alignment = {'V': best_v_align_data,
                       'D': best_d_align_data,
                       'J': best_j_align_data}

    ofile = open(task.igor_batchname+'__'+str(indexed_sequence.seq_index)+'_na.fasta', 'w')
    str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
    ofile.write(str_fasta)
    ofile.close()

    ofile = open('sequences2align.fasta', 'w')
    str_fasta = generate_str_fasta_simple(indexed_sequence, vdj_naive_alignment)
    ofile.write(str_fasta)
    ofile.close()

    import pygor3.utils as utils
    print(utils.da_heavy_pen_nuc44_vect)

    # 1. Print only sequence with corrected offset
    # 2. identify mistmatches
    # 3. show only aligment sequences

    """
    # print(seq_str)
    input_seq_align = p3.IgorAlignment_data()
    input_seq_align.seq_index = input.seq_index
    input_seq_align.strGene_seq = input.sequence
    input_seq_align.offset = 0
    input_seq_align.strGene_name = str( input.seq_index )
    input_seq_align.strGene_class = "InputSequence"
    input_seq_align.length = len(input_seq_align.strGene_seq)

    print("=" * 70)
    print(input_seq_align.strGene_seq)
    print("-" * 70)
    # print(best_v_align_data.strGene_seq)
    # print(best_d_align_data.strGene_seq)

    str_gap = '-'
    str_prefix = str_gap * (best_j_align_data.offset)
    str_tmp = str_prefix+best_j_align_data.strGene_seq
    print(str_tmp)
    print(best_j_align_data)
    # list(str_tmp)
    # from offset to offset_5_p
    str_tmp_2 = str_tmp[best_j_align_data.offset_5_p+1:best_j_align_data.offset_3_p+1]
    print(len(str_tmp_2))
    print(str_tmp_2)

    
    print("-" * 70)
    fasta_data = [input_seq_align, best_v_align_data, best_d_align_data, best_j_align_data]
    min_offset = min(fasta_data, key=lambda x: x.offset).offset
    print(min_offset)

    seq_list = list(input_seq_align.strGene_seq)
    vdj_mismatches = best_v_align_data.mismatches+best_d_align_data.mismatches+best_j_align_data.mismatches
    for ii in vdj_mismatches:
        seq_list[ii] = seq_list[ii].lower()
    input_seq_align.strGene_seq = ''.join(seq_list)
    fasta_data[0] = input_seq_align

    for item in fasta_data:
        str_prefix = 'X'*(item.offset - min_offset)
        seq_tmp = str_prefix + item.strGene_seq
        print(seq_tmp)
        aaa = seq_tmp[(item.offset_5_p - min_offset):(item.offset_3_p - min_offset)+1]
        print(item.strGene_class, item.offset_5_p, item.offset_3_p)
        print(aaa)

    # load Alignments from file


    # load Anchors with anchors extract CDR3 and naive alignment
    
    """

if __name__ == "__main__":
    main()
