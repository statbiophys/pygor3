#!/usr/bin/python3

import sqlite3

"
#-- projects table
'''
CREATE TABLE IF NOT EXISTS IgorIndexedSeq (
    seq_index integer PRIMARY KEY,
    sequence text NOT NULL,
);
'''
 
#-- tasks table
'''
CREATE TABLE IF NOT EXISTS GenomicTemplate (
    id integer PRIMARY KEY,
    name text NOT NULL,
    priority integer,
    project_id integer NOT NULL,
    status_id integer NOT NULL,
    begin_date text NOT NULL,
    end_date text NOT NULL,
    FOREIGN KEY (project_id) REFERENCES projects (id)
);
'''

def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except sqlite.Error as e:
        print(e)
 
    return conn


def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)

flnIndexedSeqCSV = "pytest/aligns/Barb02_indexed_sequences.csv"
flnIndexedSeqDB  = "test.db"#flnIndexedSeqCSV.split()

db = sqlite3.connect(flnIndexedSeqDB)

cursor = db.cursor()

# 1. Create database
# 2. Create a Table 
# 3. Read csv file line by line and insert values on database.

cursor.execute('''INSERT INTO Simulation(Directory, strF2B, strB2F, strB2B, strHrad, strBrad) VALUES(?, ?, ?, ?, ?, ?)''',
