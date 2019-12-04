#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:02:28 2019

@author: alfaceor
"""

import sqlite3
import csv
from Bio import SeqIO
import IgorAlignment_data


###################
# TODO: Generalize to other kinds of database like postgres
###################

class IgorSqliteDB:
    """
    Class to create and load table or database with sequences
    """
    def __init__(self):
        self.flnIgorSQL = "IgorDB.sql"
        self.flnIgorDB  = "" #flnIgorDB
        
        self.flnIgorIndexedSeq = ""
        
        self.flnVGeneTemplate   = ""
        self.flnDGeneTemplate   = ""
        self.flnJGeneTemplate   = ""
        
        self.flnVAlignments     = ""
        self.flnDAlignments     = ""
        self.flnJAlignments     = ""
        
        self.flnVAnchors        = ""
        self.flnJAnchors        = ""
        
        self.conn               = None
    
    #@classmethod
    def createSqliteDB(self, flnIgorDB):
        """
        Create a SQLite database with the flnIgorDB sql script.
        """
        self.flnIgorDB = flnIgorDB
        self.conn = None
        try:
            self.conn = sqlite3.connect(flnIgorDB)
            qry = open(self.flnIgorSQL, 'r').read()
            cur = self.conn.cursor()
            cur.executescript(qry)
            self.conn.commit()
            cur.close()
            #self.conn.close()
        except sqlite3.Error as e:
            print(e)
        

    def load_VDJ_Database(self, flnIgorIndexedSeq, flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, flnVAlignments, flnDAlignments, flnJAlignments):
        
        self.load_IgorIndexedSeq_FromCSV(flnIgorIndexedSeq )
        
        self.load_IgorGeneTemplate_FromFASTA("V",flnVGeneTemplate)
        self.load_IgorGeneTemplate_FromFASTA("D",flnDGeneTemplate)
        self.load_IgorGeneTemplate_FromFASTA("J",flnJGeneTemplate)
        
        self.load_IgorAlignments_FromCSV("V", flnVAlignments)
        self.load_IgorAlignments_FromCSV("D", flnDAlignments)
        self.load_IgorAlignments_FromCSV("J", flnJAlignments)
    
    ###############################################
    ####  IgorIndexedSeq Table Methods
    ###############################################
    
    def load_IgorIndexedSeq_FromCSV(self, flnIgorIndexedSeq):
        """
        Insert indexed sequence in database from csv igor indexed_seqs file.
        :param conn:
        :param csvline:
        :return:
        """
        self.flnIgorIndexedSeq = flnIgorIndexedSeq        
        
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            with open(flnIgorIndexedSeq) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    #print(csvline)
                    self.insert_IgorIndexedSeq_FromCSVline(cur, csvline)
            cur.execute('COMMIT')
            #self.conn.commit()
        except sqlite3.Error as e:
            print(e)

    
    def insert_IgorIndexedSeq_FromCSVline(self, cur, csvline):
        """
        Insert IGoR indexed_sequences on Database flnIgorDB
        :param csvline:
        """
        sql = ''' INSERT INTO IgorIndexedSeq(seq_index,sequence)
                  VALUES(?,?) '''
        
        csvline = csvline.replace('\n','')
        data = tuple(csvline.split(";"))
        if len(data) == 2:
            try:
                #cur = self.conn.cursor()
                cur.execute(sql, data)
                #self.conn.commit()
            except sqlite3.Error as e:
                print(data)
                print(e)
                pass
        

    def fetch_IgorIndexedSeq_By_seq_index(self, seq_index):
        """
        Fetch seq_index and sequence in Igor database.
        :param seq_index: string to specify the type of gene V, D or J
        :return: 
        """

        #print(gene_name)
        sqlSelect = "SELECT * FROM IgorIndexedSeq WHERE seq_index = "+str(seq_index)+";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchone()
        return (record)
    
    def fetch_IgorIndexedSeq_By_seq_indexList(self, seq_indexList):
        """
        Fetch seq_index and sequence in Igor database.
        :param seq_index: string to specify the type of gene V, D or J
        :return: 
        """

        #print(gene_name)
        strSeq_indexList = str( tuple( sorted( set( seq_indexList ) ) )  )

        sqlSelect = "SELECT * FROM IgorIndexedSeq WHERE seq_index IN "+strSeq_indexList+";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchall()
        return (record)
    
    
    
    ###############################################
    ####  IgorXGeneTemplate Tables Methods
    ###############################################
    
    def load_IgorGeneTemplate_FromFASTA(self, strGene, flnGeneTemplate):
        """
        Insert D Gene templates in database from fasta files used by IGoR.
        :param flnIgorGeneTemplate: Fasta file
        """
        # TODO: ADD A PANDAS OBJECT FOR RAPID ACCESS TO GENOMIC TEMPLATES.
        filename = {'V': self.flnVGeneTemplate, 'D' : self.flnDGeneTemplate, 'J': self.flnJGeneTemplate }
        filename[strGene] = flnGeneTemplate
        with open(filename[strGene], "r") as handle:
            for gene_id, bioRecord in enumerate(SeqIO.parse(handle, "fasta") ):
                self.insert_IgorGeneTemplate_FromBioRecord(strGene, gene_id, bioRecord)
        
   
    def insert_IgorGeneTemplate_FromBioRecord(self, strGene, gene_id, bioRecord):
        """
        Insert IGoR Gene template in Database flnIgorDB
        :param strGene: string to specify the type of gene V, D or J
        :param gene_id: id to identify the gene template
        :param bioRecord: Biopython record of the inserted sequence
        """
        sql = "INSERT INTO Igor"+strGene+"GeneTemplate("+strGene.lower()+"gene_id,gene_name,sequence) VALUES(?,?,?) "

        data = tuple([gene_id, str(bioRecord.description).strip(), str(bioRecord.seq)])
        try:
            cur = self.conn.cursor()
            cur.execute(sql, data)
            self.conn.commit()
        except sqlite3.Error as e:
            print(e)
            pass

    
    def fetch_IgorGeneTemplate_By_gene_name(self, strGene, gene_name):
        """
        Fetch Gene templates in database from fasta files used by IGoR.
        :param strGene: string to specify the type of gene V, D or J
        :param flnIgorGeneTemplate: Fasta file
        """

        sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"GeneTemplate WHERE gene_name =\""+gene_name+"\";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchone()
        return (record)        
    
    def fetch_IgorGeneTemplate_By_gene_id(self, strGene, gene_id):
        """
        Fetch Gene templates in database from fasta files used by IGoR.
        :param strGene: string to specify the type of gene V, D or J
        :param gene_id: 
        """
        sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"GeneTemplate WHERE "+strGene.lower()+"gene_id ="+str(gene_id)+";"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchone()
        return (record)        
    
    
    ###############################################
    ####  IgorXAlignments Tables Methods
    ###############################################
    
    def load_IgorAlignments_FromCSV(self, strGene, flnAlignments):
        """
        Insert Gene templates in database from fasta files used by IGoR.
        :param strGene: string to specify the type of gene V, D or J
        :param flnIgorGeneTemplate: Fasta file
        """
        strGene = strGene.upper()
        filename = {'V': self.flnVAlignments, 'D' : self.flnDAlignments, 'J': self.flnJAlignments}
        filename[strGene] = flnAlignments
        try:
            cur = self.conn.cursor()
            cur.execute('BEGIN TRANSACTION')
            with open(filename[strGene]) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    #print(csvline)
                    self.insert_IgorAlignments_FromCSVline(cur, strGene, csvline)
            
            cur.execute('COMMIT')
        except sqlite3.Error as e:
            print(e)
    
        
    
    def insert_IgorAlignments_FromCSVline(self, cur, strGene, csvline):
        """
        Insert IGoR Alignments on Database flnIgorDB
        :param strGene: string to specify the type of gene V, D or J
        :param csvline:
        """
        sqlBase = '''INSERT INTO Igor{}Alignments(
                seq_index, 
                {}gene_id, 
                score,
                offset,
                insertions,
                deletions,
                mismatches,
                length,
                offset_5_p_align,
                offset_3_p_align) 
                VALUES(?,?,?,?,?,?,?,?,?,?)
                '''
        
        strGene =  strGene.upper()
        sql = sqlBase.format(strGene.upper(), strGene.lower())
        # search on IgorGeneTemplate the corresponding Gene id 
        csvline = csvline.replace("{","[").replace("}","]").replace('\n','')
        csvlist = csvline.split(";")
        if len(csvlist) == 10 :
            gene_name = csvlist[1]
            gene_name = gene_name.strip()
            gene_id = self.fetch_IgorGeneTemplate_By_gene_name(strGene, gene_name)[0]
            #print(gene_id, gene_name)
            csvlist[1] = str(gene_id)
            data = tuple(csvlist)
            #print(data)
            try:
                #cur = self.conn.cursor()
                cur.execute(sql, data)
                #self.conn.commit()
            except sqlite3.Error as e:
                print(e)
                pass
        else:
            print(csvlist)
            
    
    def fetch_IgorAlignments_By_seq_index(self, strGene, seq_index):
        """
        Fetch IgorAlignments from database by seq_index.
        :param strGene: string to specify the type of gene V, D or J
        :param seq_index: IgorIndexedSequences index
        """
        sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"Alignments WHERE seq_index=="+str(seq_index)+" ORDER BY score DESC"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchall()
        return (record) 
    
    
    # TODO: Create an IgorAlignment_data instance with a better sql query (join the necessary tables.)
    def appendList_IgorAlignments_data_By_seq_index(self, strGene_class, seq_index, alnDataList=None):
        """
        Append to a list of IgorAlignment_data objects given gene class ("V", "D", "J"), seq_index 
        append a list to append the objects.
        :param strGene_class: string to specify the type of gene V, D or J.
        :param seq_index: IgorIndexedSequences index.
        :param alnDataList: List of IgorAlignment_data objects.
        """
        if alnDataList == None:
            alnDataList = list()
        
        try:
            alignsSQLrecords = self.fetch_IgorAlignments_By_seq_index(strGene_class, seq_index)
            for alignSQLrecord in alignsSQLrecords:
                #print(alignSQLrecord)
                gene_id = alignSQLrecord[1]
                geneTemplateRecord = self.fetch_IgorGeneTemplate_By_gene_id(strGene_class, gene_id)
                strGene_name = geneTemplateRecord[1]
                strGene_seq  = geneTemplateRecord[2]
                aln_data = IgorAlignment_data.IgorAlignment_data.load_FromSQLRecord(alignSQLrecord, strGene_name=strGene_name)
                aln_data.strGene_class = strGene_class
                aln_data.strGene_seq   = strGene_seq
                alnDataList.append(aln_data)
                #print(aln_data.strGene_name, aln_data.score, aln_data.offset, aln_data.insertions)
            return alnDataList
        except Exception as e:
            print(e)
    
    def write_IgorAlignments2Fasta(self, alnDataList):
        print(alnDataList)

        
    

    
    
    
        


#        data = tuple(csvline.split(";"))
#        print(data)
#        try:
#            cur = self.conn.cursor()
#            cur.execute(sql, data)
#            self.conn.commit()
#        except sqlite3.Error as e:
#            print(e)
#            pass
        


