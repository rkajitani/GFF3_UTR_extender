#!/usr/bin/env python

import os
import sys


class Transcript:
    def __init__(self, seq_name, strand, start, end, gene_id):
        self.seq_name = seq_name
        self.strand = strand
        self.pos_list = [(start, end)]
        self.gene_id = gene_id


class EdgeExonHit:
    def __init__(self, edge_exon_pos, ext_exon_pos, ext_trans_id):
        self.edge_exon_pos = edge_exon_pos
        self.ext_exon_pos = ext_exon_pos
        self.ext_trans_id = ext_trans_id


class ExtensionParts:
    def __init__(self):
        self.len = 0
        self.pos_list = []
        self.trans_id = ""


if len(sys.argv) != 4:
    print("usage:", os.path.basename(sys.argv[0]), "base.bed extension.bed edge_overlap.tsv", file=sys.stderr)
    sys.exit(1)


trans_dict = [dict(), dict()]
for i in (0, 1):
    # 0: base, 1: for_extension
    with open(sys.argv[i + 1]) as fin:
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            exon_id, trans_id, gene_id = f[3].split(";")
            if trans_id in trans_dict[i]:
                trans_dict[i][trans_id].pos_list.append((int(f[1]), int(f[2])))
            else:
                trans_dict[i][trans_id] = Transcript(f[0], f[5], int(f[1]), int(f[2]), gene_id)
    for trans_info in trans_dict[i].values():
        trans_info.pos_list.sort(key=lambda x: x[0])


hit_dict = [dict(), dict()]
with open(sys.argv[3]) as fin:
    for ln in fin:
        f = ln.rstrip("\n").split("\t")
        edge_exon_pos = (int(f[1]), int(f[2]))
        ext_exon_pos = (int(f[7]), int(f[8]))
        ovl_len = int(f[12])
        if edge_exon_pos[1] - edge_exon_pos[0] != ovl_len:
            continue
        edge_exon_id, edge_trans_id, edge_gene_id = f[3].split(";")
        ext_exon_id, ext_trans_id, ext_gene_id = f[9].split(";")

        direction = 0
        # 0: left, 1: right
        if edge_exon_id.split(".")[-1] == "right":
            direction = 1

        oposite_dir = direction ^ 1
        if edge_exon_pos[oposite_dir] == ext_exon_pos[oposite_dir]:
            if edge_trans_id in hit_dict:
                hit_dict[direction][edge_trans_id].append(EdgeExonHit(edge_exon_pos, ext_exon_pos, ext_trans_id))
            else:
                hit_dict[direction][edge_trans_id] = [EdgeExonHit(edge_exon_pos, ext_exon_pos, ext_trans_id)]


edge_max_ext_parts = dict()
for direction in (0, 1):
    oposite_dir = direction ^ 1
    for edge_trans_id, hit_list in hit_dict[direction].items():
        for hit in hit_list:
            if edge_trans_id not in edge_max_ext_parts:
                edge_max_ext_parts[edge_trans_id] = [ExtensionParts(), ExtensionParts()]
            max_ext_parts = edge_max_ext_parts[edge_trans_id][direction]

            ext_exon_pos_list = trans_dict[1][hit.ext_trans_id].pos_list
            ext_len = 0
            pos_list = []
            if direction == 0:
                for i in range(0, len(ext_exon_pos_list)):
                    ext_exon_pos = ext_exon_pos_list[i]
                    if hit.edge_exon_pos[oposite_dir] == ext_exon_pos[oposite_dir]:
                        ext_len += hit.edge_exon_pos[direction] - ext_exon_pos[direction]
                        ext_len += sum([x[1] - x[0] for x in ext_exon_pos_list[:i]])
                        pos_list += ext_exon_pos_list[:i]
                        if ext_exon_pos[direction] < hit.edge_exon_pos[direction]:
                            pos_list.append((ext_exon_pos[direction], hit.edge_exon_pos[direction]))
                        break
            else:
                for i in reversed(range(0, len(ext_exon_pos_list))):
                    ext_exon_pos = ext_exon_pos_list[i]
                    if hit.edge_exon_pos[oposite_dir] == ext_exon_pos[oposite_dir]:
                        ext_len += ext_exon_pos[direction] - hit.edge_exon_pos[direction]
                        ext_len += sum([x[1] - x[0] for x in ext_exon_pos_list[i + 1:]])
                        if hit.edge_exon_pos[direction] < ext_exon_pos[direction]:
                            pos_list.append((hit.edge_exon_pos[direction], ext_exon_pos[direction]))
                        pos_list += ext_exon_pos_list[i + 1:]
                        break

            if max_ext_parts.len < ext_len:
                max_ext_parts.len = ext_len
                max_ext_parts.pos_list = pos_list.copy()
                max_ext_parts.trans_id = hit.ext_trans_id


for center_trans_id, center_trans_info in trans_dict[0].items():
    print(center_trans_info.seq_name, center_trans_info.gene_id, center_trans_id, center_trans_info.strand, sep="\t", end="\t")
    print(sum([x[1] - x[0] for x in center_trans_info.pos_list]), end="\t")
    print(";".join([str(x[0]) for x in center_trans_info.pos_list]), end="\t")
    print(";".join([str(x[1]) for x in center_trans_info.pos_list]), end="")
    for i in (0, 1):
        if center_trans_id in edge_max_ext_parts and edge_max_ext_parts[center_trans_id][i].len > 0:
            ext_parts = edge_max_ext_parts[center_trans_id][i]
            print("", ext_parts.trans_id, trans_dict[1][ext_parts.trans_id].strand, sep="\t", end="\t")
            print(sum([x[1] - x[0] for x in ext_parts.pos_list]), end="\t")
            print(";".join([str(x[0]) for x in ext_parts.pos_list]), end="\t")
            print(";".join([str(x[1]) for x in ext_parts.pos_list]), end="")
        else:
            print("\t", "\t".join(["."] * 5), sep="", end="")
    print()
