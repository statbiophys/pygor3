#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:18:45 2019

@author: alfaceor
"""

class IgorAlignment_data:
    def __init__(self):
        self.seq_index     = -1
        self.gene_id       = -1
        self.score         = -1
        self.offset        = 0
        self.insertions    = list()
        self.deletions     = list() 
        self.mismatches    = list()
        self.length        = 0
        self.offset_5_p    = 0
        self.offset_3_p    = 0
        
        self.strGene_name  = ""
        self.strGene_class = ""
        self.strGene_seq   = ""

    def to_dict(self):
        dictAlignment_data       =  {
            "seq_index"     : self.seq_index  , \
            "gene_id"       : self.gene_id    , \
            "score"         : self.score      , \
            "offset"        : self.offset     , \
            "insertions"    : self.insertions , \
            "deletions"     : self.deletions  , \
            "mismatches"    : self.mismatches , \
            "length"        : self.length     , \
            "offset_5_p"    : self.offset_5_p , \
            "offset_3_p"    : self.offset_3_p , \
            "strGene_name"  : self.strGene_name    , \
            "strGene_class" : self.strGene_class   , \
            "strGene_seq"   : self.strGene_seq
            }
        
        return dictAlignment_data
    
#    def __str__(self):
#        strIgorAlignment_data = "[" + \
#            str(self.seq_index)  + ", " + \
#            str(self.gene_id)    + ", " + \
#            str(self.score)      + ", " + \
#            str(self.offset)     + ", " + \
#            str(self.insertions) + ", " + \
#            str(self.deletions)  + ", " + \
#            str(self.mismatches) + ", " + \
#            str(self.length)     + ", " + \
#            str(self.offset_5_p) + ", " + \
#            str(self.offset_3_p) + ", " + \
#            self.strGene_name    + ", " + \
#            self.strGene_class   + ", " + \
#            self.strGene_seq     + ", " + \
#            "]"
#            
#        #return strIgorAlignment_data
#        return str( self.to_dict() )
    
    @classmethod
    def load_FromSQLRecord(cls, sqlRecordAlign, strGene_name=""):
        """
        Return a IgorAlignment_data instance from a IgorSqlRecord.
        :param sqlRecordAlign: record of a sql database table.
        :param strGene_name: gene_name associated to the record.
        :return: IgorAlignment_data instance
        """
        cls = IgorAlignment_data()
        try:
            cls.seq_index    = int  (sqlRecordAlign[0])
            cls.gene_id      = int  (sqlRecordAlign[1])
            cls.score        = float(sqlRecordAlign[2])
            cls.offset       = int  (sqlRecordAlign[3])
            cls.insertions   = eval (sqlRecordAlign[4])
            cls.deletions    = eval (sqlRecordAlign[5])
            cls.mismatches   = eval (sqlRecordAlign[6])
            cls.length       = int  (sqlRecordAlign[7])
            cls.offset_5_p   = int  (sqlRecordAlign[8])
            cls.offset_3_p   = int  (sqlRecordAlign[9])
            # TODO: Bestway to retrieve the name of the gene_name
            if strGene_name == None:
                cls.strGene_name = str(cls.gene_id)
            else:
                cls.strGene_name = strGene_name
            return cls
        except Exception as e:
            print(e)
            raise e
        
#   
#    @classmethod
#    def writeFastaAlignments(sqlRecordIgorIndexedSeq, sqlRecordAlignList):
#        """
#        
#        """
#        seq_index = sqlRecordIgorIndexedSeq[0]
#        strSeq    = sqlRecordIgorIndexedSeq[1]
#        print("> "+str(seq_index))
#        print(strSeq)
#
#        
#def fromManySqlRecords2AlignmentList(sqlRecordAlignments):
#    if len(sqlRecordAlignments)>1:
#        igorAlndata = IgorAlignment_data()
#        for sqlRecordAlign in sqlRecordAlignments:
#            if len(sqlRecordAlignments) == 10:
#                igorAlndata.load_FromSQLRecord(sqlRecordAlign)
#
