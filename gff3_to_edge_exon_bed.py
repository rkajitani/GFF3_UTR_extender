#!/usr/bin/env python

import sys
import gff_utils


class Transcript:
    def __init__(self, seq_name, strand, start, end, gene_id):
        self.seq_name = seq_name
        self.strand = strand
        self.pos_list = [(start, end)]
        self.gene_id = gene_id


if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "in.gff3 exon_feature_str", file=sys.stderr)
    sys.exit(1)

exon_feat_str = sys.argv[2]

gff3_tree = gff_utils.GFF3Tree()
gff3_tree.set_relations(sys.argv[1])

transcript_dict = dict()
with open(sys.argv[1]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        feat_type = f[2]
        if feat_type != exon_feat_str:
            continue
        seq_name = f[0]
        start = int(f[3])
        end = int(f[4])
        strand = f[6]
        attr_str = f[8]
        feat_id = gff_utils.gff3_attr_get_id(attr_str)
        parent_id = gff3_tree.get_parent_id(feat_id)
        root_id = gff3_tree.get_root_id(feat_id)
        if parent_id in transcript_dict:
            transcript_dict[parent_id].pos_list.append((start, end))
        else:
            transcript_dict[parent_id] = Transcript(seq_name, strand, start, end, root_id)

for trans_id, transcript in transcript_dict.items():
    sorted_pos = sorted(transcript.pos_list, key=lambda x: x[0])
    print(transcript.seq_name, sorted_pos[0][0] - 1, sorted_pos[0][1], f"{trans_id}.left;{trans_id};{transcript.gene_id}", 0, transcript.strand, sep="\t")
    print(transcript.seq_name, sorted_pos[-1][0] - 1, sorted_pos[-1][1], f"{trans_id}.right;{trans_id};{transcript.gene_id}", 0, transcript.strand, sep="\t")
