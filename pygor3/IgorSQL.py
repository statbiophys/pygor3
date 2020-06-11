# SQL commands to create database
sqlcmd_ct = dict()

sqlcmd_ct['indexed_sequences'] = """
-- IgorIndexedSeq table
CREATE TABLE IF NOT EXISTS IgorIndexedSeq (
    seq_index integer PRIMARY KEY,
    sequence text NOT NULL
);
"""

sqlcmd_ct['indexed_CDR3'] = """
-- IgorIndexedSeq table
CREATE TABLE IF NOT EXISTS IgorIndexedCDR3 (
    seq_index integer PRIMARY KEY,
    v_anchor integer,
    j_anchor integer
);
"""
# """
# -- IgorIndexedSeq table
# CREATE TABLE IF NOT EXISTS IgorIndexedCDR3 (
#     seq_index integer PRIMARY KEY,
#     v_anchor integer,
#     j_anchor integer,
#     CDR3nt text,
#     CDR3aa text
# );
# """

sqlcmd_ct['genomicVs'] = """
-- IgorVGeneTemplate table
CREATE TABLE IF NOT EXISTS IgorVGeneTemplate (
    vgene_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);
"""

sqlcmd_ct['genomicDs'] = """
-- IgorDGeneTemplate table
CREATE TABLE IF NOT EXISTS IgorDGeneTemplate (
    dgene_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);
"""

sqlcmd_ct['genomicJs'] = """
-- IgorJGeneTemplate table
CREATE TABLE IF NOT EXISTS IgorJGeneTemplate (
    jgene_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);
"""

sqlcmd_ct['geneVCDR3Anchors'] = """
-- IgorVGeneCDR3Anchors table
CREATE TABLE IF NOT EXISTS IgorVGeneCDR3Anchors (
    vgene_id integer PRIMARY KEY,
    gene_name text,
    anchor_index integer,
    function text
);
"""

sqlcmd_ct['geneJCDR3Anchors'] = """
-- IgorJGeneCDR3Anchors table
CREATE TABLE IF NOT EXISTS IgorJGeneCDR3Anchors (
    jgene_id integer PRIMARY KEY,
    gene_name text,
    anchor_index integer,
    function text
);
"""

sqlcmd_ct['V_alignments'] = """
-- IgorVAlignments table
CREATE TABLE IF NOT EXISTS IgorVAlignments (
    seq_index integer NOT NULL,
    vgene_id integer NOT NULL,
    score real,
    offset integer,
    insertions text NOT NULL,
    deletions  text NOT NULL,
    mismatches text NOT NULL,
    length integer,
    offset_5_p_align integer,
    offset_3_p_align integer,
    FOREIGN KEY (seq_index) REFERENCES IgorIndexedSeq   (seq_index),
    FOREIGN KEY (vgene_id)  REFERENCES IgorVGeneTemplate(vgene_id)
    --PRIMARY KEY (seq_index, vgene_id)
);
"""

sqlcmd_ct['D_alignments'] = """
-- IgorDAlignments table
CREATE TABLE IF NOT EXISTS IgorDAlignments (
    seq_index integer NOT NULL,
    dgene_id integer NOT NULL,
    score real,
    offset integer,
    insertions text NOT NULL,
    deletions  text NOT NULL,
    mismatches text NOT NULL,
    length integer,
    offset_5_p_align integer,
    offset_3_p_align integer,
    FOREIGN KEY (seq_index) REFERENCES IgorIndexedSeq   (seq_index),
    FOREIGN KEY (dgene_id)  REFERENCES IgorDGeneTemplate(dgene_id)
    --PRIMARY KEY (seq_index, dgene_id)
);
"""

sqlcmd_ct['J_alignments'] = """
-- IgorJAlignments table
CREATE TABLE IF NOT EXISTS IgorJAlignments (
    seq_index integer NOT NULL,
    jgene_id integer NOT NULL,
    score real,
    offset integer,
    insertions text NOT NULL,
    deletions  text NOT NULL,
    mismatches text NOT NULL,
    length integer,
    offset_5_p_align integer,
    offset_3_p_align integer,
    FOREIGN KEY (seq_index) REFERENCES IgorIndexedSeq   (seq_index),
    FOREIGN KEY (jgene_id)  REFERENCES IgorJGeneTemplate(jgene_id)
    --PRIMARY KEY (seq_index, jgene_id)
);
"""


# XXX: New tables use these.
sqlcmd_ct['XXXXXXXXXXXXXXXX'] = """
YYYYYYYYYYYYYYYYYYYYYYY
"""

sqlcmd_ct['XXXXXXXXXXXXXXXX'] = """
YYYYYYYYYYYYYYYYYYYYYYY
"""

sqlcmd_ct['XXXXXXXXXXXXXXXX'] = """
YYYYYYYYYYYYYYYYYYYYYYY
"""

