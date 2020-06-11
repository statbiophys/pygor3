-- IgorIndexedSeq table
CREATE TABLE IF NOT EXISTS IgorIndexedSeq (
    seq_index integer PRIMARY KEY,
    sequence text NOT NULL
);


-- IgorVGeneTemplate table
CREATE TABLE IF NOT EXISTS IgorVGeneTemplate (
    vgene_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);

-- IgorDGeneTemplate table
CREATE TABLE IF NOT EXISTS IgorDGeneTemplate (
    dgene_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);

-- IgorJGeneTemplate table
CREATE TABLE IF NOT EXISTS IgorJGeneTemplate (
    jgene_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);

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

-- IgorDAlignments table
CREATE TABLE IF NOT EXISTS IgorDAlignments (
    seq_index integer NOT NULL,
    dgene_id integer NOT NULL,
    score real,
    offset integer,
    insertions text NOT NULL,
    deletions  text NOT NULL,
    mismatches text NOT NULL,
    n_insertions integer,
    n_deletions integer,
    n_mismatches integer,
    length integer,
    offset_5_p_align integer,
    offset_3_p_align integer,
    FOREIGN KEY (seq_index) REFERENCES IgorIndexedSeq   (seq_index),
    FOREIGN KEY (dgene_id)  REFERENCES IgorDGeneTemplate(dgene_id)
    --PRIMARY KEY (seq_index, dgene_id)
);

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

