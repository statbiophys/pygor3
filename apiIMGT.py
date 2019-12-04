#!/usr/bin/python3
from __future__ import print_function
import os
import IgorModel

import requests
from bs4 import BeautifulSoup
from Bio import SeqIO

try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3


# extract sequence from the IMGT website
# \return raw fasta data file.
def getIMGT_fasta(url):
    """ Extract fasta information from the typical IMGT format
    """
    print(url)
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    result = soup.find_all('pre')[1].text
    return result


# \param seqDescritption : is the header description in fasta file
# \return key : a generated key to create a dictionary
# \return startPos : the start position of the sequence based on ACCesion number sequence.
def genKey(seqDescription):
    fields = seqDescription.split("|")
    key    = "|".join(fields[:4])
    return key

# return the start position of sequence based on the fasta file description
def getStartPos(seqDescription):
    fields = seqDescription.split("|")
    #startPosLIGMDB = map(int, fields[5].split(".."))[0]
    startPosLIGMDB = int ( fields[5].split("..")[0] )
    return startPosLIGMDB

# return a list of Sequence records.
def genSeqRecords(url):
    rawfasta     = getIMGT_fasta(url)
    return list(SeqIO.parse(StringIO(rawfasta), "fasta")) # FIXME: Is better to preserve the generator but ... you know maybe later XD

# return a dictinary like key = fastaStartPos
def genAnchDict(url):
    recordsAnch = genSeqRecords(url)
    seqsDictAnch = {}
    for rec in recordsAnch:
        seqsDictAnch[genKey(rec.description)] =  getStartPos(rec.description)
    #print("genAnchDict", seqsDictAnch)
    return seqsDictAnch


def write_V_TemplateFiles(gene, specie,):
    DIR_REF_GENOME = specie+"/"+gene+"/ref_genome/"
    gene = gene + "V"
    print ("DIR_REF_GENOME = ", DIR_REF_GENOME)
    # 1. Get the Alleles/Gene sequences from the IMGT website (IMGT/GENE-DB)
    URL_IMGT_GENEDB = "http://www.imgt.org/genedb/GENElect?"
    QUERY_GENE      = "query=7.2+"+gene+"&species="+specie
    QUERY_V_2CYS    = "query=8.1+"+gene+"&species="+specie+"&IMGTlabel=2nd-CYS"
    urlVgene        = URL_IMGT_GENEDB + QUERY_GENE
    urlV_2CYS       = URL_IMGT_GENEDB + QUERY_V_2CYS
    
    # 2. load data from IMGT url 
    records = genSeqRecords(urlVgene)
    records = map(lambda x: x.upper(), records)
    
    # 3. Write it in a fasta file.
    flnGenomicVs = DIR_REF_GENOME+"genomicVs.fasta"
    with open(flnGenomicVs, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")
    
    # generate dictionaries with the anchors
    dictV_2CYS = genAnchDict(urlV_2CYS)
    #print( len(dictV_2CYS), len(records) )
    # write the files 
    flnAnchors = DIR_REF_GENOME+"V_gene_CDR3_anchors.csv"
    ofileAnch  = open(flnAnchors, "w")
    ofileAnch.write("gene;anchor_index"+"\n")
    CSVDELIM = ";"
    for rec in records:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posV_2CYS = -1
        if key in dictV_2CYS.keys():
            posV_2CYS =  dictV_2CYS[key] - getStartPos(rec.description) 
            posAnch  = posV_2CYS
            ofileAnch.write(rec.description+CSVDELIM+str(posAnch)+"\n")
    
        else:
            print("No anchor is found for : "+ rec.description)
    
    ofileAnch.close()


def write_D_TemplateFiles(gene, specie):
    DIR_REF_GENOME = specie+"/"+gene+"/ref_genome/"
    gene = gene + "D"
    # 1. Get the Alleles/Gene sequences from the IMGT website (IMGT/GENE-DB)
    URL_IMGT_GENEDB = "http://www.imgt.org/genedb/GENElect?"
    QUERY_GENE      = "query=7.2+"+gene+"&species="+specie
    urlDgene        = URL_IMGT_GENEDB + QUERY_GENE
    
    # 2. load data from IMGT url 
    records = genSeqRecords(urlDgene)
    records = map(lambda x: x.upper(), records)
    
    # 3. Write it in a fasta file.
    flnGenomicDs = DIR_REF_GENOME+"genomicDs.fasta"
    with open(flnGenomicDs, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")



def write_J_TemplateFiles(gene, specie):
#    gene    = "TRBJ"
#    species = "Homo+sapiens"
    DIR_REF_GENOME = specie+"/"+gene+"/ref_genome/"
    gene = gene + "J"
    # 1. Get the Alleles/Gene sequences from the IMGT website (IMGT/GENE-DB)
    URL_IMGT_GENEDB = "http://www.imgt.org/genedb/GENElect?"
    QUERY_GENE      = "query=7.2+"+gene+"&species="+specie
    QUERY_J_PHE     = "query=8.1+"+gene+"&species="+specie+"&IMGTlabel=J-PHE"
    QUERY_J_TRP     = "query=8.1+"+gene+"&species="+specie+"&IMGTlabel=J-TRP"
    urlJgene        = URL_IMGT_GENEDB + QUERY_GENE
    urlJ_PHE        = URL_IMGT_GENEDB + QUERY_J_PHE
    urlJ_TRP        = URL_IMGT_GENEDB + QUERY_J_TRP
    
    # 2. load data from IMGT url 
    records = genSeqRecords(urlJgene)
    records = map(lambda x: x.upper(), records)
    
    # 3. Write it in a fasta file.
    flnGenomicJs = DIR_REF_GENOME+"genomicJs.fasta"
    with open(flnGenomicJs, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")
    
    # generate dictionaries with the anchors
    dictJ_PHE = genAnchDict(urlJ_PHE)
    dictJ_TRP = genAnchDict(urlJ_TRP)
    # write the files 
    flnAnchors = DIR_REF_GENOME+"J_gene_CDR3_anchors.csv"
    ofileAnch  = open(flnAnchors, "w")
    ofileAnch.write("gene;anchor_index"+"\n")
    CSVDELIM = ";"
    for rec in records:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posJ_PHE = -1
        posJ_TRP = -1
        if key in dictJ_PHE.keys():
            posJ_PHE =  dictJ_PHE[key] - getStartPos(rec.description) 
            posAnch  = posJ_PHE
            ofileAnch.write(rec.description+CSVDELIM+str(posAnch)+"\n")
    
        # FIXME: this should be elif, but I want to check If the system is consistent
        elif key in dictJ_TRP.keys():
            posJ_TRP = dictJ_TRP[key] -  getStartPos(rec.description) 
            posAnch  = posJ_TRP
            ofileAnch.write(rec.description+CSVDELIM+str(posAnch)+"\n")
        else:
            print("No anchor is found for : "+ rec.description)
    
    ofileAnch.close()
    
    



#gene    = "TRBJ"
#specie = "Homo+sapiens"
#write_J_TemplateFiles("TRBJ", specie)
#write_V_TemplateFiles("TRBV", specie)
#write_D_TemplateFiles("TRBD", specie)

def writeModelVDJ_Parms(gene, species, model_parms_file=None):
    DIR_REF_GENOME = species+"/"+gene+"/ref_genome/"
    DIR_MODELS     = species+"/"+gene+"/models/"
    flnGenomicVs = DIR_REF_GENOME+"genomicVs.fasta"
    flnGenomicDs = DIR_REF_GENOME+"genomicDs.fasta"
    flnGenomicJs = DIR_REF_GENOME+"genomicJs.fasta"
    Vgenes = dict()
    Dgenes = dict()
    Jgenes = dict()
    
    if model_parms_file == None :
        model_parms_file = DIR_MODELS +"model_parms.txt"
    else:
        model_parms_file=model_parms_file
    
    for record in list(SeqIO.parse(flnGenomicVs, "fasta") ):
        Vgenes[record.description] = record.seq
    
    for record in list(SeqIO.parse(flnGenomicDs, "fasta") ):
        Dgenes[record.description] = record.seq
    
    for record in list(SeqIO.parse(flnGenomicJs, "fasta") ):
        Jgenes[record.description] = record.seq
    
    with open(DIR_MODELS +"model_parms.txt", "w") as fw:
        print("@Event_list", file=fw)
        strPriority="7"
        print("#GeneChoice;V_gene;Undefined_side;"+strPriority+";v_choice", file=fw)
        for iv, v in enumerate(Vgenes):
            print("%{};{};{}".format(v, Vgenes[v], iv), file=fw)

        strPriority="6"
        print("#GeneChoice;D_gene;Undefined_side;"+strPriority+";d_gene", file=fw)
        for iD, D in enumerate(Dgenes):
            print("%{};{};{}".format(D, Dgenes[D], iD), file=fw)

        strPriority="7"
        print("#GeneChoice;J_gene;Undefined_side;"+strPriority+";j_choice", file=fw)
        for ij, j in enumerate(Jgenes):
            print("%{};{};{}".format(j, Jgenes[j], ij), file=fw)

        # for deletion/insertion/..., default values
        strPriority="5"
        print("#Deletion;V_gene;Three_prime;"+strPriority+";v_3_del", file=fw)
        delVmin=-4
        delVmax=17
        for iDelV, DelV in enumerate(range(delVmin, delVmax)):
            print("%{};{}".format(DelV, iDelV), file=fw)

        strPriority="5"
        print("#Deletion;D_gene;Three_prime;"+strPriority+";d_3_del", file=fw)

        delD3min=-4
        delD3max=17
        for iDelD, DelD in enumerate(range(delD3min, delD3max)):
            print("%{};{}".format(DelD, iDelD), file=fw)

        strPriority="5"
        print("#Deletion;D_gene;Five_prime;"+strPriority+";d_5_del", file=fw)
        delD5min=-4
        delD5max=17
        for iDelD, DelD in enumerate(range(delD5min, delD5max)):
            print("%{};{}".format(DelD, iDelD), file=fw)

        strPriority="5"
        print("#Deletion;J_gene;Five_prime;"+strPriority+";j_5_del", file=fw)
        delJmin=-4
        delJmax=19
        for iDelJ, DelJ in enumerate(range(delJmin, delJmax)):
            print("%{};{}".format(DelJ, iDelJ), file=fw)

        strPriority="4"
        print("#Insertion;VD_genes;Undefined_side;"+strPriority+";vd_ins", file=fw)
        insVD=31
        for i in range(insVD):
            print("%{0};{0}".format(i), file=fw)
         
        strPriority="2"
        print("#Insertion;DJ_gene;Undefined_side;"+strPriority+";dj_ins", file=fw)
        insDJ=31
        for i in range(insDJ):
            print("%{0};{0}".format(i), file=fw)

        strPriority="3"
        print("#DinucMarkov;VD_genes;Undefined_side;"+strPriority+";vd_dinucl", file=fw)
        for iN, N in enumerate(["A", "C", "G", "T"]):
            print("%{};{}".format(N, iN), file=fw)
        strPriority="1"
        print("#DinucMarkov;DJ_gene;Undefined_side;"+strPriority+";dj_dinucl", file=fw)
        for iN, N in enumerate(["A", "C", "G", "T"]):
            print("%{};{}".format(N, iN), file=fw)

        # edges, default again
        print(("@Edges\n"
               "%GeneChoice_V_gene_Undefined_side_prio7_size{0};"
               "Deletion_V_gene_Three_prime_prio5_size21\n"
               "%GeneChoice_D_gene_Undefined_side_prio6_size{1};"
               "Deletion_D_gene_Three_prime_prio5_size21\n"
               "%GeneChoice_D_gene_Undefined_side_prio6_size{1};"
               "Deletion_D_gene_Five_prime_prio5_size21\n"
               "%GeneChoice_J_gene_Undefined_side_prio7_size{2};"
               "Deletion_J_gene_Five_prime_prio5_size23\n"
               "%GeneChoice_J_gene_Undefined_side_prio7_size{2};"
               "GeneChoice_D_gene_Undefined_side_prio6_size{1}\n"
               "%Deletion_D_gene_Five_prime_prio5_size21;"
               "Deletion_D_gene_Three_prime_prio5_size21\n"
               "@ErrorRate\n#SingleErrorRate").format(
                   len(Vgenes), len(Dgenes), len(Jgenes)), file=fw)

        error_rate = 0.001  # default value
        print(error_rate, file=fw)
        
#    parms = IgorModel.Model_Parms(model_parms_file=DIR_MODELS +"model_parms.txt")
#    parms.genMarginalFile(model_marginals_file=DIR_MODELS +"model_marginals.txt")


def writeModelVDJ_Marginals(gene, species, model_parms_file=None, model_marginals_file=None):
    DIR_MODELS     = species+"/"+gene+"/models/"
    if model_parms_file == None :
        model_parms_file = DIR_MODELS +"model_parms.txt"
        print(model_parms_file)
        parms = IgorModel.Model_Parms(model_parms_file=model_parms_file)
    else:
        parms = IgorModel.Model_Parms(model_parms_file=model_parms_file)
    
    if model_marginals_file == None :
        model_marginals_file = DIR_MODELS +"model_marginals.txt"
    else:
        model_marginals_file = model_marginals_file
    
    with open(model_marginals_file, "w") as fw:
        # v_choice
        strEvent = "v_choice"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+"]\n")
        fw.write("#\n")
        fw.write("%")
        unifProb = (1./Dim_v)
        for vv in range(Dim_v):
            fw.write(str(unifProb))
            if not (vv == Dim_v-1):
                fw.write(",")
        fw.write("\n")
        
        # j_choice
        strEvent="j_choice"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict['v_choice'])
        Dim_j = len(parms.Event_dict['j_choice'])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_j)+"]\n")
        unifProb = 1./(Dim_j)
        for vv in range(Dim_v):
            fw.write("#[v_choice,"+str(vv)+"]\n")
            fw.write("%")
            for jj in range(Dim_j):
                fw.write(str(unifProb))
                if not (jj == Dim_j-1):
                    fw.write(",")
            fw.write("\n")
        
        #d_gene
        strEvent = "d_gene"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict['v_choice'])
        Dim_j = len(parms.Event_dict['j_choice'])
        Dim_d = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_j)+","+str(Dim_d)+"]\n")
        unifProb = 1./(Dim_d)
        for vv in range(Dim_v):
            for jj in range(Dim_j):
                fw.write("#[v_choice,"+str(vv)+"],[j_choice,"+str(jj)+"]\n")
                fw.write("%")
                for dd in range(Dim_d):
                    fw.write(str(unifProb))
                    if not (dd == Dim_d-1):
                        fw.write(",")
                fw.write("\n")
        
        #v_3_del
        strEvent = "v_3_del"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict['v_choice'])
        Dim_v3del = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_v3del)+"]\n")
        unifProb = 1./(Dim_v3del)
        for vv in range(Dim_v):
            fw.write("#[v_choice,"+str(vv)+"]\n")
            fw.write("%")
            for v3 in range(Dim_v3del):
                fw.write(str(unifProb))
                if not (v3 == Dim_v3del-1):
                    fw.write(",")
            fw.write("\n")
        
        #d_5_del
        strEvent = "d_5_del"
        fw.write("@"+strEvent+"\n")
        Dim_d = len(parms.Event_dict['d_gene'])
        Dim_d5del = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_d)+","+str(Dim_d5del)+"]\n")
        unifProb = 1./(Dim_d5del)
        for dd in range(Dim_d):
            fw.write("#[d_gene,"+str(dd)+"]\n")
            fw.write("%")
            for d5 in range(Dim_d5del):
                fw.write(str(unifProb))
                if not (d5 == Dim_d5del-1):
                    fw.write(",")
            fw.write("\n")
        
        #j_5_del
        strEvent = "j_5_del"
        fw.write("@"+strEvent+"\n")
        Dim_j = len(parms.Event_dict['j_choice'])
        Dim_j5del = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_j5del)+"]\n")
        unifProb = 1./(Dim_j5del)
        for jj in range(Dim_j):
            fw.write("#[j_choice,"+str(jj)+"]\n")
            fw.write("%")
            for j5 in range(Dim_j5del):
                fw.write(str(unifProb))
                if not (j5 == Dim_j5del-1):
                    fw.write(",")
            fw.write("\n")
        
        #vd_ins
        strEvent = "vd_ins"
        fw.write("@"+strEvent+"\n")
        Dim_vdins = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_vdins)+"]\n")
        unifProb = 1./(Dim_vdins)
        fw.write("#\n")
        fw.write("%")
        for vd in range(Dim_vdins):
            fw.write(str(unifProb))
            if not (vd == Dim_vdins-1):
                fw.write(",")
        fw.write("\n")
        
        #vd_dinucl
        strEvent = "vd_dinucl"
        fw.write("@"+strEvent+"\n")
        Dim_vddinu = len(parms.Event_dict[strEvent])*len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_vddinu)+"]\n")
        unifProb = 1./(Dim_vddinu)
        fw.write("#\n")
        fw.write("%")
        for vn in range(Dim_vddinu):
            fw.write(str(unifProb))
            if not (vn == Dim_vddinu-1):
                fw.write(",")
        fw.write("\n")
        
        #dj_ins
        strEvent = "dj_ins"
        fw.write("@"+strEvent+"\n")
        Dim_djins = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_djins)+"]\n")
        unifProb = 1./(Dim_djins)
        fw.write("#\n")
        fw.write("%")
        for dj in range(Dim_djins):
            fw.write(str(unifProb))
            if not (dj == Dim_djins-1):
                fw.write(",")
        fw.write("\n")
        
        #dj_dinucl
        strEvent = "dj_dinucl"
        fw.write("@"+strEvent+"\n")
        Dim_djdinu = len(parms.Event_dict[strEvent])*len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_djdinu)+"]\n")
        unifProb = 1./(Dim_djdinu)
        fw.write("#\n")
        fw.write("%")
        for dn in range(Dim_djdinu):
            fw.write(str(unifProb))
            if not (dn == Dim_djdinu-1):
                fw.write(",")
        fw.write("\n")


def writeModelVJ_Parms(gene, species, model_parms_file=None):
    DIR_REF_GENOME = species+"/"+gene+"/ref_genome/"
    DIR_MODELS     = species+"/"+gene+"/models/"
    flnGenomicVs = DIR_REF_GENOME+"genomicVs.fasta"
    flnGenomicJs = DIR_REF_GENOME+"genomicJs.fasta"
    Vgenes = dict()
    Jgenes = dict()
    
    if model_parms_file == None :
        model_parms_file = DIR_MODELS +"model_parms.txt"
        parms = IgorModel.Model_Parms(model_parms_file=model_parms_file)
    else:
        parms = IgorModel.Model_Parms(model_parms_file=model_parms_file)
    
    for record in list(SeqIO.parse(flnGenomicVs, "fasta") ):
        Vgenes[record.description] = record.seq
    
    for record in list(SeqIO.parse(flnGenomicJs, "fasta") ):
        Jgenes[record.description] = record.seq
    
    with open(model_parms_file, "w") as fw:
        print("@Event_list", file=fw)
        print("#GeneChoice;V_gene;Undefined_side;7;v_choice", file=fw)
        for iv, v in enumerate(Vgenes):
            print("%{};{};{}".format(v, Vgenes[v], iv), file=fw)
        print("#GeneChoice;J_gene;Undefined_side;7;j_choice", file=fw)
        for ij, j in enumerate(Jgenes):
            print("%{};{};{}".format(j, Jgenes[j], ij), file=fw)

        # for deletion/insertion/..., default values
        print("#Deletion;V_gene;Three_prime;5;v_3_del", file=fw)
        for iDelV, DelV in enumerate(range(-4, 17)):
            print("%{};{}".format(DelV, iDelV), file=fw)
        print("#Deletion;J_gene;Five_prime;5;j_5_del", file=fw)
        for iDelJ, DelJ in enumerate(range(-4, 19)):
            print("%{};{}".format(DelJ, iDelJ), file=fw)
        print("#Insertion;VJ_genes;Undefined_side;4;vj_ins", file=fw)
        for i in range(31):
            print("%{0};{0}".format(i), file=fw)
        print("#DinucMarkov;VJ_genes;Undefined_side;3;vj_dinucl", file=fw)
        for iN, N in enumerate(["A", "C", "G", "T"]):
            print("%{};{}".format(N, iN), file=fw)

        # edges, default again
        print(("@Edges\n"
               "%GeneChoice_V_gene_Undefined_side_prio7_size{0};"
               "GeneChoice_J_gene_Undefined_side_prio7_size{1};\n"
               "%GeneChoice_V_gene_Undefined_side_prio7_size{0};"
               "Deletion_V_gene_Three_prime_prio5_size21\n"
               "%GeneChoice_J_gene_Undefined_side_prio7_size{1};"
               "Deletion_J_gene_Five_prime_prio5_size23\n"
               "@ErrorRate\n#SingleErrorRate").format(
                   len(Vgenes), len(Jgenes)), file=fw)

        error_rate = 0.001  # default value
        print(error_rate, file=fw)
    
    
#    parms = IgorModel.Model_Parms(model_parms_file=DIR_MODELS +"model_parms.txt")
#    parms.genMarginalFile(model_marginals_file=DIR_MODELS +"model_marginals.txt")


def writeModelVJ_Marginals(gene, species, model_parms_file=None, model_marginals_file=None):
    DIR_MODELS     = species+"/"+gene+"/models/"
    if model_parms_file == None :
        model_parms_file = DIR_MODELS +"model_parms.txt"
        print(model_parms_file)
        parms = IgorModel.Model_Parms(model_parms_file=model_parms_file)
    else:
        parms = IgorModel.Model_Parms(model_parms_file=model_parms_file)
    
    if model_marginals_file == None :
        model_marginals_file = DIR_MODELS +"model_marginals.txt"
    else:
        model_marginals_file = model_marginals_file
    
    with open(model_marginals_file, "w") as fw:
        # v_choice
        strEvent = "v_choice"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+"]\n")
        fw.write("#\n")
        fw.write("%")
        unifProb = (1./Dim_v)
        for vv in range(Dim_v):
            fw.write(str(unifProb))
            if not (vv == Dim_v-1):
                fw.write(",")
        fw.write("\n")
        
        # j_choice
        strEvent="j_choice"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict['v_choice'])
        Dim_j = len(parms.Event_dict['j_choice'])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_j)+"]\n")
        unifProb = 1./(Dim_j)
        for vv in range(Dim_v):
            fw.write("#[v_choice,"+str(vv)+"]\n")
            fw.write("%")
            for jj in range(Dim_j):
                fw.write(str(unifProb))
                if not (jj == Dim_j-1):
                    fw.write(",")
            fw.write("\n")
        
        #v_3_del
        strEvent = "v_3_del"
        fw.write("@"+strEvent+"\n")
        Dim_v = len(parms.Event_dict['v_choice'])
        Dim_v3del = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_v3del)+"]\n")
        unifProb = 1./(Dim_v3del)
        for vv in range(Dim_v):
            fw.write("#[v_choice,"+str(vv)+"]\n")
            fw.write("%")
            for v3 in range(Dim_v3del):
                fw.write(str(unifProb))
                if not (v3 == Dim_v3del-1):
                    fw.write(",")
            fw.write("\n")
        
        #j_5_del
        strEvent = "j_5_del"
        fw.write("@"+strEvent+"\n")
        Dim_j = len(parms.Event_dict['j_choice'])
        Dim_j5del = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_v)+","+str(Dim_j5del)+"]\n")
        unifProb = 1./(Dim_j5del)
        for jj in range(Dim_j):
            fw.write("#[j_choice,"+str(jj)+"]\n")
            fw.write("%")
            for j5 in range(Dim_j5del):
                fw.write(str(unifProb))
                if not (j5 == Dim_j5del-1):
                    fw.write(",")
            fw.write("\n")
        
        #vj_ins
        strEvent = "vj_ins"
        fw.write("@"+strEvent+"\n")
        Dim_vjins = len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_vjins)+"]\n")
        unifProb = 1./(Dim_vjins)
        fw.write("#\n")
        fw.write("%")
        for vj in range(Dim_vjins):
            fw.write(str(unifProb))
            if not (vj == Dim_vjins-1):
                fw.write(",")
        fw.write("\n")
        
        #vj_dinucl
        strEvent = "vj_dinucl"
        fw.write("@"+strEvent+"\n")
        Dim_vjdinu = len(parms.Event_dict[strEvent])*len(parms.Event_dict[strEvent])
        fw.write("$Dim["+str(Dim_vjdinu)+"]\n")
        unifProb = 1./(Dim_vjdinu)
        fw.write("#\n")
        fw.write("%")
        for vjn in range(Dim_vjdinu):
            fw.write(str(unifProb))
            if not (vjn == Dim_vjdinu-1):
                fw.write(",")
        fw.write("\n")
        
    

def makeDirectories(gene, species):
    os.system("mkdir -p "+species)
    os.system("mkdir -p "+species+"/"+gene)
    os.system("mkdir -p "+species+"/"+gene+"/ref_genome")
    os.system("mkdir -p "+species+"/"+gene+"/ref_genome")
    os.system("mkdir -p "+species+"/"+gene+"/models")
    


