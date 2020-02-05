

flnIndexedSeqCSV = "pytest/aligns/Barb02_indexed_sequences.csv"
flnIndexedSeqDB  = "test.db"#flnIndexedSeqCSV.split()

#db = sqlite3.connect(flnIndexedSeqDB)
#cursor = db.cursor()

# 1. Create database
conn = create_connection(flnIndexedSeqDB)
# 2. Create a Table 
create_table(conn, create_table_sql)

#csvline="seq_index;sequence"
#csvline="11;CGCACACAGCAGGAGGACTCCGCCGTGTATCTCTGTGCCAGCAGCTTACGAGGGGGAGGCAGCAATCAGCCCCAGCATTTTGGTGAT"
#csvline="7;CGCACACAGCAGGAGGACTCCGCCGTGTATCTCTGTGCCAGCAGCTTACGAGGGGGAGGCAGCAATCAGCCCCAGCATTTTGGTGAT"
# 3. Read csv file line by line and insert values on database.
with open(flnIndexedSeqCSV) as fp:
    csvline = fp.readline()
    while csvline:
#        try:
        print(csvline)
        insert_IgorIndexedSeq_FromCSVline(conn, csvline)
        csvline = fp.readline()
#        except sqlite3.Error as e:
#            print(e)
#            pass
            

conn.close()
#print("database closed")

