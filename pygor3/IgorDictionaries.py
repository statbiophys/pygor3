import copy
# FIXME: PUT THIS IN __init__.py ? is a good practice?
# option wrappers for IGoR
# IGoR command line options
igor_option_path_dict={
    "alpha": "tcr_alpha",
    "beta" : "tcr_beta",
    "light": "light",
    "heavy_naive" : "bcr_heavy",
    "heavy_memory": "bcr_heavy"
}

# Options used to align in  IGoR
tmp_dict_options = {
    '---thresh': {'active': False, 'value': '15', 'dict_options': {}},
    '---matrix': {'active': False, 'value': 'path/to/file', 'dict_options': {}},
    '---gap_penalty': {'active': False, 'value': 'X', 'dict_options': {}},
    '---best_align_only': {'active': False, 'value': '', 'dict_options': {}}
}
igor_align_dict_options = {
    '--all':
        {'active': True, 'value': '',
         'dict_options':
             copy.deepcopy(tmp_dict_options)
         },
    '--V':
        {'active': False, 'value': '',
         'dict_options':
             copy.deepcopy(tmp_dict_options)
         },
    '--D':
        {'active': False, 'value': '',
         'dict_options':
             copy.deepcopy(tmp_dict_options)
         },
    '--J':
        {'active': False, 'value': '',
         'dict_options':
             copy.deepcopy(tmp_dict_options)
         }
}
igor_align_dict_options['--V']['dict_options']['---thresh']['value'] = '50'


igor_evaluate_dict_options = {
    '--N_iter':
        {'active': False, 'value': '5'},
    '--L_thresh':
        {'active': False, 'value': ''},
    '--P_ratio_thresh':
        {'active': False, 'value': '0.0'}, # Check all possible scenarios
    '--MLSO':
        {'active': False, 'value': ''},
    '--infer_only':
        {'active': False, 'value': ''},
    '--not_infer':
        {'active': False, 'value': ''},
    '--fix_err':
        {'active': False, 'value': ''},

}


igor_infer_dict_options = {
    '--N_iter':
        {'active': False, 'value': '5'},
    '--L_thresh':
        {'active': False, 'value': ''},
    '--P_ratio_thresh':
        {'active': False, 'value': '0.0'}, # Check all possible scenarios
    '--MLSO':
        {'active': False, 'value': ''},
    '--infer_only':
        {'active': False, 'value': ''},
    '--not_infer':
        {'active': False, 'value': ''},
    '--fix_err':
        {'active': False, 'value': ''},

}
# --N_iter N
# --L_thresh X
# --P_ratio_thresh X
# --MLSO
# --infer_only eventnickname1 eventnickname2
# --not_infer eventnickname1
# --fix_err

# Options used to output in  IGoR
igor_output_dict_options = {
    '--scenarios':
        {'active': True, 'value': '10',
         'dict_options':
             {}
         },
    '--Pgen':
        {'active': True, 'value': '',
         'dict_options':
             {}
         },
    '--coverage':
        {'active': False, 'value': '',
         'dict_options':
             {}
         }
}

igor_generate_dict_options = {
    '--noerr':
        {'active': False, 'value': None,
         'dict_options': {}
        },
    '--CDR3':
        {'active': False, 'value': None,
         'dict_options': {}
        },
    '--name':
        {'active': False, 'value': None,
         'dict_options': {}
        },
    '--seed':
        {'active': False, 'value': None,
         'dict_options': {}
        }
}

# FIXME: MAKE A DICTIONARY
igor_file_id_list = [
# aligns directory
'indexed_sequences',
'indexed_CDR3',

'aligns_V_alignments',
'aligns_D_alignments',
'aligns_J_alignments',

# inference directory
'infer_final_parms',
'infer_final_marginals',

# evaluate directory
'evaluate_final_parms',
'evaluate_final_marginals',

# output directory
'output_pgen',
'output_scenarios',
'output_coverage',

# IGoR ref_genome
'genomicVs',
'genomicDs',
'genomicJs',

'V_gene_CDR3_anchors',
'J_gene_CDR3_anchors',

# IGoR models
'model_parms',
'model_marginals',

]

# dictionary for the status of the filenames used in IGoR at different stages.
igor_batch_dict = {
    # aligns directory
    'indexed_sequences' :  {'filename': "", 'status': False},
    'indexed_CDR3' :  {'filename': "", 'status': False},

    'aligns_V_alignments' :  {'filename': "", 'status': False},
    'aligns_D_alignments' :  {'filename': "", 'status': False},
    'aligns_J_alignments' :  {'filename': "", 'status': False},

    # inference directory
    'infer_final_parms' :  {'filename': "", 'status': False},
    'infer_final_marginals' :  {'filename': "", 'status': False},

    # evaluate directory
    'evaluate_final_parms' :  {'filename': "", 'status': False},
    'evaluate_final_marginals' :  {'filename': "", 'status': False},

    # output directory
    'output_pgen' :  {'filename': "", 'status': False},
    'output_scenarios' :  {'filename': "", 'status': False},
    'output_coverage' :  {'filename': "", 'status': False},

    # IGoR ref_genome
    'genomicVs' :  {'filename': "", 'status': False},
    'genomicDs' :  {'filename': "", 'status': False},
    'genomicJs' :  {'filename': "", 'status': False},

    'V_gene_CDR3_anchors' :  {'filename': "", 'status': False},
    'J_gene_CDR3_anchors' :  {'filename': "", 'status': False},

    # IGoR models
    'model_parms' :  {'filename': "", 'status': False},
    'model_marginals' :  {'filename': "", 'status': False}

}

def update_igor_batch_dict(igor_wd, igor_batchname):
    tmp_prefix_aligns = igor_wd + "/aligns/" + igor_batchname
    igor_batch_dict['indexed_sequences']['filename'] = tmp_prefix_aligns + "_indexed_sequences.csv"
    igor_batch_dict['indexed_CDR3']['filename'] = tmp_prefix_aligns + "_indexed_CDR3.csv"
    igor_batch_dict['aligns_V_alignments']['filename'] = tmp_prefix_aligns + "_V_alignments.csv"
    igor_batch_dict['aligns_D_alignments']['filename'] = tmp_prefix_aligns + "_D_alignments.csv"
    igor_batch_dict['aligns_J_alignments']['filename'] = tmp_prefix_aligns + "_J_alignments.csv"

    tmp_prefix = igor_wd + "/" + igor_batchname
    igor_batch_dict['infer_final_parms']['filename'] = tmp_prefix + "_inference/" + "final_parms.txt"
    igor_batch_dict['infer_final_marginals']['filename'] = tmp_prefix + "_inference/" + "final_marginals.txt"
    igor_batch_dict['evaluate_final_parms']['filename'] = tmp_prefix + "_evaluate/" + "final_parms.txt"
    igor_batch_dict['evaluate_final_marginals']['filename'] = tmp_prefix + "_evaluate/" + "final_marginals.txt"
    igor_batch_dict['output_pgen']['filename'] = tmp_prefix + "_output/" + "Pgen_counts.csv"
    igor_batch_dict['output_scenarios']['filename'] = tmp_prefix + "_output/" + "best_scenarios_counts.csv"
    igor_batch_dict['output_coverage']['filename'] = tmp_prefix + "_output/" + "coverage.csv"

igor_options = {}
