#!/usr/bin/env python3
import pygor3 as p3


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
    db.load_VDJ_Database(
        task.igor_fln_indexed_sequences,
        genomes.fln_genomicVs,
        genomes.fln_genomicDs,
        genomes.fln_genomicJs,
        task.igor_fln_align_V_alignments,
        task.igor_fln_align_D_alignments,
        task.igor_fln_align_J_alignments)

    seq_index = 1
    seq_str = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
    print("seq_str", seq_str)

    alnDataListV = db.fetch_IgorAlignments_By_seq_index('V', seq_index)
    alnDataListD = db.fetch_IgorAlignments_By_seq_index('D', seq_index)
    alnDataListJ = db.fetch_IgorAlignments_By_seq_index('J', seq_index)

    # Only best alignments
    v_align_data = p3.IgorAlignment_data.load_FromSQLRecord(alnDataListV[0])
    d_align_data = p3.IgorAlignment_data.load_FromSQLRecord(alnDataListD[0])
    j_align_data = p3.IgorAlignment_data.load_FromSQLRecord(alnDataListJ[0])

    best_v_align_data = db.fetch_best_IgorAlignments_By_seq_index('V', seq_index)



    print("v_align_data")
    print(v_align_data.to_dict())
    print("align_data")

    seq_offset = 0
    min_offset = min([seq_offset, v_align_data.offset, d_align_data.offset, j_align_data.offset])
    if min_offset < 0:
        str_seq_offset = "-"*(-min_offset)
        base_offset = -min_offset
    else:
        str_seq_offset = ""
        base_offset = seq_offset

    v_gene_seq = db.fetch_IgorGeneTemplate_By_gene_id('V', v_align_data.gene_id)[2]
    d_gene_seq = db.fetch_IgorGeneTemplate_By_gene_id('D', d_align_data.gene_id)[2]
    j_gene_seq = db.fetch_IgorGeneTemplate_By_gene_id('J', j_align_data.gene_id)[2]

    print('='*80)
    # print("v_align_data", v_align_data)
    v_align_data.strGene_class = 'V'
    gene_id, v_align_data.strGene_name, v_align_data.strGene_seq = db.fetch_IgorGeneTemplate_By_gene_id(v_align_data.strGene_class, v_align_data.gene_id)
    # print("v_align_data", v_align_data)

    print(v_gene_seq)
    print(v_gene_seq[:v_align_data.offset_3_p])
    aaa = "X"*(-v_align_data.offset)
    print(aaa+seq_str)
    print(aaa + seq_str[: v_align_data.offset_3_p])
    print(aaa + seq_str[: v_align_data.offset_3_p])

    # genomes.dict_genomicVs[]
    # v_align_data.strGene_name =

    bbb = db.fetch_IgorGeneTemplate_By_gene_id('V', v_align_data.gene_id)
    print(v_align_data.gene_id, bbb)
    v_align_data.strGene_class = 'V'
    v_align_data.strGene_name = bbb[1]
    v_align_data.strGene_seq = bbb[2]
    print(v_align_data.to_dict())

    # print(alnDataListV)


    # print(genomes_dict)

    # print(genomes.df_genomicVs[ genomes.df_genomicVs['name'] == align_data.strGene_name])

    # print(type(genomes.df_genomicVs[genomes.df_genomicVs['name'] == align_data.strGene_name]['value'].values ))

    # FIXME : USING JUST FILES WITHOUT DATABASE
    # open and read alignment file line by line
    ofile = open(task.igor_fln_align_V_alignments, 'r')
    header = ofile.readline()
    print(header)
    line = ofile.readline()
    line = ofile.readline()
    line = ofile.readline()
    print(line)

    # line = ofile.readline()
    # print(line)
    align_data = p3.IgorAlignment_data.load_FromCSVLine(line)
    print(align_data.to_dict())
    print(align_data.strGene_name)

    genomes_dict = (genomes.df_genomicVs.set_index('name').to_dict())['value']
    print(type(genomes_dict[align_data.strGene_name]))
    print(align_data.offset_5_p)
    print(align_data.offset_3_p)
    print(len(align_data.mismatches))

    print(align_data.to_dict())
    print(line)

    ofile.close()


    best_v_align_data2 = db.get_best_IgorAlignment_data_By_seq_index('V', seq_index)
    best_d_align_data2 = db.get_best_IgorAlignment_data_By_seq_index('D', seq_index)
    best_j_align_data2 = db.get_best_IgorAlignment_data_By_seq_index('J', seq_index)
    naive_vdj_alignment = [ ]

    print("="*70)
    print(seq_str)
    input_seq_align = p3.IgorAlignment_data()
    input_seq_align.seq_index = seq_index
    input_seq_align.strGene_seq = seq_str
    input_seq_align.offset = 0
    input_seq_align.strGene_name = str( seq_index )
    input_seq_align.strGene_class = "InputSequence"
    input_seq_align.length = len(input_seq_align.strGene_seq)

    print(input_seq_align)
    print(best_v_align_data2)
    print(best_d_align_data2)
    print(best_j_align_data2)

    print("-" * 70)
    fasta_data = [input_seq_align, best_v_align_data2, best_d_align_data2, best_j_align_data2]
    min_offset = min(fasta_data, key=lambda x: x.offset).offset
    print(min_offset)

    seq_list = list(input_seq_align.strGene_seq)
    vdj_mismatches = best_v_align_data2.mismatches+best_d_align_data2.mismatches+best_j_align_data2.mismatches
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

if __name__ == "__main__":
    main()
