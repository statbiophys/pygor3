/*
-- IgorIndexedSeq table
CREATE TABLE IF NOT EXISTS IgorIndexedSeq (
    seq_index integer PRIMARY KEY,
    sequence text NOT NULL
);
*/


-- IgorBestScenarios table
CREATE TABLE IF NOT EXISTS IgorDBBestScenariosVDJ (
    seq_index integer NOT NULL,
    scenario_rank integer,
    scenario_proba_cond_seq real,
    id_v_choice integer,
    id_j_choice integer,
    id_d_gene integer,
    id_v_3_del integer,
    id_d_5_del integer,
    id_d_3_del integer,
    id_j_5_del integer,
    id_vd_ins integer,
    vd_dinucl text,
    id_dj_ins integer,
    dj_dinucl text,
    mismatches text,
    mismatcheslen integer
);


-- MP_v_choice table
CREATE TABLE IF NOT EXISTS MP_v_choice (
    event_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);

-- MP_j_choice table
CREATE TABLE IF NOT EXISTS MP_j_choice (
    event_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);

-- MP_d_gene table
CREATE TABLE IF NOT EXISTS MP_d_gene (
    event_id integer PRIMARY KEY,
    gene_name text NOT NULL,
    sequence text NOT NULL
);

-- MP_v_3_del table
CREATE TABLE IF NOT EXISTS MP_v_3_del (
    event_id integer PRIMARY KEY,
    event_value text NOT NULL
);


-- MP_d_5_del table
CREATE TABLE IF NOT EXISTS MP_d_5_del (
    event_id integer PRIMARY KEY,
    event_value text NOT NULL
);


-- MP_d_3_del table
CREATE TABLE IF NOT EXISTS MP_d_3_del (
    event_id integer PRIMARY KEY,
    event_value text NOT NULL
);

-- MP_j_5_del table
CREATE TABLE IF NOT EXISTS MP_j_5_del (
    event_id integer PRIMARY KEY,
    event_value text NOT NULL
);

-- MP_vd_ins table
CREATE TABLE IF NOT EXISTS MP_vd_ins (
    event_id integer PRIMARY KEY,
    event_value text NOT NULL
);

-- MP_vd_dinucl table
CREATE TABLE IF NOT EXISTS MP_vd_dinucl (
    event_id integer PRIMARY KEY,
    event_value text NOT NULL
);


