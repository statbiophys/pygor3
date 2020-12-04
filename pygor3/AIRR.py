
"""
J gene with allele.
If referring to a known reference sequence
in a database the relevant gene/allele nomenclature
should be followed (e.g., IGHJ4*02 if using IMGT/GENE-DB).
"""
import collections

class AIRR_gene:
    def __init__(self, gene_type):
        self.gene_type = gene_type # v, d, d2, j
        self.call = None
        self.score = None
        self.identity = None
        self.support = None
        self.cigar = None

        self.sequence_start = None
        self.sequence_end = None
        self.germline_start = None
        self.germline_end = None
        self.alignment_start = None
        self.alignment_end = None

    def to_dict(self):
        dicto = collections.OrderedDict()
        ptype = self.gene_type+"_"
        dicto[ptype + "call"] = self.call
        dicto[ptype + "score"] = self.score
        dicto[ptype + "identity"] = self.identity
        dicto[ptype + "support"] = self.support
        dicto[ptype + "cigar"] = self.cigar

        dicto[ptype + "sequence_start"] = self.sequence_start
        dicto[ptype + "sequence_end"] = self.sequence_end
        dicto[ptype + "germline_start"] = self.germline_start
        dicto[ptype + "germline_end"] = self.germline_end
        dicto[ptype + "alignment_start"] = self.alignment_start
        dicto[ptype + "alignment_end"] = self.alignment_end

        return dicto


class AIRR_VDJ_rearrangement:
    def __init__(self, sequence_id=None, sequence=None, scenario_rank=None, scenario_proba_cond_seq=None, pgen=None):
        self.sequence_id = sequence_id
        self.sequence = sequence
        self.scenario_rank = scenario_rank
        self.scenario_proba_cond_seq = scenario_proba_cond_seq
        self.pgen = pgen

        self.quality = None
        self.sequence_aa = None
        self.rev_comp = None
        self.productive = None
        self.vj_in_frame = None
        self.stop_codon = None
        self.complete_vdj = None
        self.locus = None

        #genes data
        self.v_data = AIRR_gene("v")
        self.d_data = AIRR_gene("d")
        self.j_data = AIRR_gene("j")


        self.sequence_alignment = None
        self.quality_alignment = None
        self.sequence_alignment_aa = None

        self.junction = None
        self.junction_aa = None

        self.n1_length = None

        self.np1 = None
        self.np1_aa = None
        self.np1_length = None

        self.n2_length = None

        self.np2 = None
        self.np2_aa = None
        self.np2_length = None

        self.p3v_length = None
        self.p5d_length = None
        self.p3d_length = None
        self.p5j_length = None

    def to_dict(self):
        dicto = collections.OrderedDict()
        dicto["sequence_id"] = self.sequence_id
        dicto["sequence"] = self.sequence
        dicto["scenario_rank"] = self.scenario_rank
        dicto["scenario_proba_cond_seq"] = self.scenario_proba_cond_seq
        dicto["pgen"] = self.pgen

        dicto.update(self.v_data.to_dict())
        dicto.update(self.d_data.to_dict())
        dicto.update(self.j_data.to_dict())

        dicto["quality"] = self.quality
        dicto["sequence_aa"] = self.sequence_aa
        dicto["rev_comp"] = self.rev_comp
        dicto["productive"] = self.productive
        dicto["vj_in_frame"] = self.vj_in_frame
        dicto["stop_codon"] = self.stop_codon
        dicto["complete_vdj"] = self.complete_vdj
        dicto["locus"] = self.locus
        dicto["sequence_alignment"] = self.sequence_alignment
        dicto["quality_alignment"] = self.quality_alignment
        dicto["sequence_alignment_aa"] = self.sequence_alignment_aa
        dicto["junction"] = self.junction
        dicto["junction_aa"] = self.junction_aa

        dicto["n1_length"] = self.n1_length

        dicto["np1"] = self.np1
        dicto["np1_aa"] = self.np1_aa
        dicto["np1_length"] = self.np1_length

        dicto["n2_length"] = self.n2_length

        dicto["np2"] = self.np2
        dicto["np2_aa"] = self.np2_aa
        dicto["np2_length"] = self.np2_length

        dicto["p3v_length"] = self.p3v_length
        dicto["p5d_length"] = self.p5d_length
        dicto["p3d_length"] = self.p3d_length
        dicto["p5j_length"] = self.p5j_length

        return dicto

    @staticmethod
    def list_of_fields():
        cls = AIRR_VDJ_rearrangement()
        return list(cls.to_dict().keys())



