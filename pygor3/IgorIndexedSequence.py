#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 12:24:34 2019

@author: alfaceor
"""
import pandas as pd

class IgorIndexedSequence:
    def __init__(self):
        self.seq_index     = -1
        self.sequence      = ""
        
    def to_dict(self):
        dictIndexedSequence       =  {
            "seq_index"     : self.seq_index  , \
            "sequence"      : self.sequence
            }
        
        return dictIndexedSequence
    
    @classmethod
    def load_FromCSVline(cls, csvline, delimiter=";"):
        cls = IgorIndexedSequence()
        csvsplit = csvline.replace("\n", "").split(";")
        try:
            cls.seq_index    = int  (csvsplit[0])
            cls.sequence     = int  (csvsplit[1])
        except Exception as e:
            print(e)
            raise e
        return cls
    
    @classmethod
    def load_FromSQLRecord(cls, sqlRecord):
        cls = IgorIndexedSequence()
        try:
            cls.seq_index    = int  (sqlRecord[0])
            cls.sequence     = int  (sqlRecord[1])
        except Exception as e:
            print(e)
            raise e
        return cls


