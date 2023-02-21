#!/usr/bin/env python

import sys
import gff_utils


if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "in.gff3", file=sys.stderr)
    sys.exit(1)

gff3_tree = gff_utils.GFF3Tree()
gff3_tree.set_relations(sys.argv[1])

transcript_dict = dict()
with open(sys.argv[1]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat_type = f[2]
        if feat_type != "exon":
            continue
        seq_name = f[0]
        start = int(f[3])
        end = int(f[4])
        strand = f[6]
        attr_str = f[8]
        feat_id = gff_utils.gff3_attr_get_id(attr_str)
        parent_id = gff3_tree.get_parent_id(feat_id)
        print(seq_name, start - 1, end, f"{feat_id};{parent_id}", 0, strand, sep="\t")
