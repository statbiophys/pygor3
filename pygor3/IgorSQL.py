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
    j_anchor integer,
    CDR3 text,
    CDR3_aa text
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
    vgene_id integer,
    anchor_index integer,
    function text,
    FOREIGN KEY (vgene_id)  REFERENCES IgorVGeneTemplate(vgene_id)
);
"""

sqlcmd_ct['geneJCDR3Anchors'] = """
-- IgorJGeneCDR3Anchors table
CREATE TABLE IF NOT EXISTS IgorJGeneCDR3Anchors (
    jgene_id integer,
    anchor_index integer,
    function text,
    FOREIGN KEY (jgene_id)  REFERENCES IgorJGeneTemplate(jgene_id)
);
"""

sqlcmd_ct['V_alignments'] = """
-- IgorVAlignments table
CREATE TABLE IF NOT EXISTS IgorVAlignments (
    seq_index integer NOT NULL,
    vgene_id integer NOT NULL,
    score integer,
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
    score integer,
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
    score integer,
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


############# IGOR MODEL TABLES #############
sqlcmd_ct['MP_Event_list'] = """
-- IgorMP_Event_list table
--     event_id integer,
CREATE TABLE IF NOT EXISTS IgorMP_Event_list (
    nickname text NOT NULL PRIMARY KEY,
    event_type text,
    seq_type text,
    seq_side text,
    priority integer,
    realizations_table text,
    name text
);
"""

sqlcmd_ct['ER_event_template'] = """
-- ER_event_template table
CREATE TABLE IF NOT EXISTS {} (
    id integer NOT NULL PRIMARY KEY,
    value text,
    name text
);
"""

sqlcmd_ct['MP_Edges'] = """
-- IgorMP_Edges table
CREATE TABLE IF NOT EXISTS IgorMP_Edges (
    parent_event text NOT NULL,
    child_event text NOT NULL,
    FOREIGN KEY (parent_event) REFERENCES IgorMP_Event_list (nickname),
    FOREIGN KEY (child_event)  REFERENCES IgorMP_Event_list (nickname)
);
"""

sqlcmd_ct['MP_ErrorRate'] = """
-- IgorMP_ErrorRate table
CREATE TABLE IF NOT EXISTS IgorMP_ErrorRate (
    error_type text NOT NULL,
    error_values text NOT NULL
);
"""

def sqlcmd_ct_Model_Marginals(event_nickname, list_dependencies:list):
    # from the xarray get the list of events something like:
    lista = [event_nickname]+list_dependencies
    str_column_ct = "id_{} integer NOT NULL" #.format(event_nickname)
    str_foreign_key_ct = "FOREIGN KEY (id_{}) REFERENCES IgorER_{} (id)" #.format(event_nickname)

    sqlcmd_table_fields_ct = ",\n".join([str_column_ct.format(evento_nickname) for evento_nickname in lista])
    sqlcmd_foreign_keys_ct = ",\n".join([str_foreign_key_ct.format(evento_nickname, evento_nickname) for evento_nickname in lista])

    sqlcmd_ct_aux = """
            -- MM_XXXXXX table
            CREATE TABLE IF NOT EXISTS IgorMM_{} (
                -- Events id columns
                {},
                P real,
                -- Foreign keys
                {}
            );
            """

    sqlcmd_ct = sqlcmd_ct_aux.format(event_nickname, sqlcmd_table_fields_ct, sqlcmd_foreign_keys_ct)
    return sqlcmd_ct

def sqlcmd_ct_Model_Marginals_DinucMarkov(event_nickname, list_dependencies:list):
    # from the xarray get the list of events something like:
    # lista = [event_nickname]+list_dependencies
    lista = list_dependencies
    str_column_ct = "id_{} integer NOT NULL" #.format(event_nickname)
    str_foreign_key_ct = "FOREIGN KEY (id_{}) REFERENCES IgorER_{} (id)" #.format(event_nickname)

    sqlcmd_table_fields_ct = ",\n".join([str_column_ct.format(evento_nickname) for evento_nickname in lista])
    sqlcmd_foreign_keys_ct = ",\n".join([str_foreign_key_ct.format(evento_nickname, event_nickname) for evento_nickname in lista])

    sqlcmd_ct_aux = """
            -- MM_XXXXXX table
            CREATE TABLE IF NOT EXISTS IgorMM_{} (
                -- Events id columns
                {},
                P real,
                -- Foreign keys
                {}
            );
            """

    sqlcmd_ct = sqlcmd_ct_aux.format(event_nickname, sqlcmd_table_fields_ct, sqlcmd_foreign_keys_ct)
    return sqlcmd_ct

############## IGOR BEST SCENARIOS ##############
sqlcmd_ct['Pgen'] = """
-- IgorPgen table
CREATE TABLE IF NOT EXISTS IgorPgen (
    seq_index integer PRIMARY KEY,
    Pgen_estimate real NOT NULL
);
"""


sqlcmd_ct['BestScenario_template'] = """
-- IgorBestScenarios table
CREATE TABLE IF NOT EXISTS IgorBestScenarios (
    seq_index integer NOT NULL,
    scenario_rank integer NOT NULL,
    scenario_proba_cond_seq real NOT NULL,
    {},
    mismatches text,
    mismatcheslen integer,
    {}
);
"""

def sqlcmd_ct_BestScenarios(nickname_event_type_list:list):
    """
    param nickname_event_type_list: list of tuples (nickname, event_type)
    return sql command to create BestScenarios table.
    """
    # FIXME: ID EVENTS EXCEPT WITH DinucMarkov
    # SELECT * from IgorMP_Event_list; or SELECT nickname, event_type from IgorMP_Event_list;
    # SELECT COUNT(*) FROM IgorMP_Edges WHERE child_event =='d_gene';
    # SELECT COUNT(IgorMP_Edges.parent_event) AS nparents, IgorMP_Event_list.nickname FROM IgorMP_Edges, IgorMP_Event_list WHERE IgorMP_Edges.child_event = IgorMP_Event_list.nickname GROUP BY IgorMP_Event_list.nickname;

    str_column_ct_with_id = "id_{} integer NOT NULL"
    str_column_ct_DinucMarkov = "{} text"
    str_foreign_key_ct = "FOREIGN KEY (id_{}) REFERENCES IgorER_{} (id)"

    str_column_list = list()
    str_foreign_key_list = list()
    for nickname, event_type in nickname_event_type_list:
        if event_type == 'DinucMarkov':
            str_column_list.append(str_column_ct_DinucMarkov.format(nickname))
        else:
            str_column_list.append(str_column_ct_with_id.format(nickname))
            str_foreign_key_list.append(str_foreign_key_ct.format(nickname, nickname))

    str_columns = ",\n".join(str_column_list)
    str_foreign_keys = ",\n".join(str_foreign_key_list)

    sqlcmd_ct_aux = """
    -- Best_scenarios table
    CREATE TABLE IF NOT EXISTS IgorBestScenarios (
        seq_index integer NOT NULL,
        scenario_rank integer NOT NULL,
        scenario_proba_cond_seq real NOT NULL,
        {},
        mismatches text,
        mismatcheslen integer,
        {}
    );
    """

    return sqlcmd_ct_aux.format(str_columns, str_foreign_keys)

# usage : sql_cmd_template_qry_cols.format(tablename='IgorMM_d_3_del')
sqlcmd_template_qry_cols = "SELECT name FROM pragma_table_info('{tablename}');"

# sqlcmd_qry = dict()
# sqlcmd_qry['IgorMarginals'] =

# sqlcmd_ct['ER_XXXXXX'] = """
# -- MP_Event_list table
# CREATE TABLE IF NOT EXISTS MM_XXXXXX (
#     XXXXX_id integer NOT NULL,
#     realization_value text,
#     realization_name text
# );
# """



# # XXX: New tables use these.
# sqlcmd_ct['XXXXXXXXXXXXXXXX'] = """
# YYYYYYYYYYYYYYYYYYYYYYY
# """
#
# sqlcmd_ct['XXXXXXXXXXXXXXXX'] = """
# YYYYYYYYYYYYYYYYYYYYYYY
# """
#
# sqlcmd_ct['XXXXXXXXXXXXXXXX'] = """
# YYYYYYYYYYYYYYYYYYYYYYY
# """
#
# sqlcmd_ct['MM_XXXXXX'] = """
# -- MM_XXXXXX table
# CREATE TABLE IF NOT EXISTS MM_XXXXXX (
#     XXXXX_id integer NOT NULL,
#     P real,
#     FOREIGN KEY (XXXXX_id) REFERENCES ER_XXXXXX (id),
# );
# """

# --------------------------- SQL TABLENAME PATTERNS --------------------------- #

sql_tablename_patterns = ['IgorIndexedSeq',
                          'IgorIndexedCDR3', 'Igor%GeneTemplate', 'Igor%GeneCDR3Anchors', 'Igor%Alignments',
                          'IgorMP_%', 'IgorER_%', 'IgorMM_%',
                          'IgorPgen', 'IgorBestScenarios']

sql_tablename_patterns_dict = dict()
sql_tablename_patterns_read_seqs = ['IgorIndexedSeq']
sql_tablename_patterns_ref_genome = ['Igor%GeneTemplate', 'Igor%GeneCDR3Anchors']
sql_tablename_patterns_align = ['IgorIndexedCDR3', 'Igor%Alignments']
sql_tablename_patterns_model = ['IgorMP_%', 'IgorER_%', 'IgorMM_%']
sql_tablename_patterns_output = ['IgorPgen', 'IgorBestScenarios']
sql_tablename_patterns_pgen = ['IgorPgen']
sql_tablename_patterns_scenarios = ['IgorBestScenarios']
sql_tablename_patterns_dict['read_seqs'] = sql_tablename_patterns_read_seqs
sql_tablename_patterns_dict['ref_genome'] = sql_tablename_patterns_ref_genome
sql_tablename_patterns_dict['align'] = sql_tablename_patterns_align
sql_tablename_patterns_dict['model'] = sql_tablename_patterns_model
sql_tablename_patterns_dict['output'] = sql_tablename_patterns_output
sql_tablename_patterns_dict['pgen'] = sql_tablename_patterns_pgen
sql_tablename_patterns_dict['scenarios'] = sql_tablename_patterns_scenarios

# ----------------------------------------------------------------------------- #


