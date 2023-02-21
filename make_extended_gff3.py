#!/usr/bin/env python

import os
import sys


if len(sys.argv) != 2:
    print("usage:", os.path.basename(sys.argv[0]), "extension_info.tsv", file=sys.stderr)
    sys.exit(1)


source_name = "lemon"

print("##gff-version 3")
with open(sys.argv[1]) as fin:
    for ln in fin:
        f = ln.rstrip("\n").split("\t")
        seq_name = f[0]
        gene_id = f[1]
        trans_id = f[2]
        strand = f[3]

        if f[5] != ".":
            center_pos_list = [(int(x[0]), int(x[1])) for x in zip(f[5].split(";"), f[6].split(";"))]
        else:
            center_pos_list = []

        if f[10] != ".":
            left_pos_list = [(int(x[0]), int(x[1])) for x in zip(f[10].split(";"), f[11].split(";"))]
        else:
            left_pos_list = []

        if f[15] != ".":
            right_pos_list = [(int(x[0]), int(x[1])) for x in zip(f[15].split(";"), f[16].split(";"))]
        else:
            right_pos_list = []

        center_pos_list.sort(key=lambda x: x[0])
        left_pos_list.sort(key=lambda x: x[0])
        right_pos_list.sort(key=lambda x: x[0])

        raw_list = left_pos_list + center_pos_list + right_pos_list
        exon_pos_list = []
        pre_pos = raw_list[0]
        for i in range(1, len(raw_list)):
            if pre_pos[1] >= raw_list[i][0]:
                pre_pos = (pre_pos[0], raw_list[i][1])
            else:
                exon_pos_list.append(pre_pos)
                pre_pos = raw_list[i]
        exon_pos_list.append(pre_pos)

        print(
            seq_name,
            source_name,
            "gene",
            exon_pos_list[0][0] + 1,
            exon_pos_list[-1][1],
            ".",
            strand,
            ".",
            f"ID={gene_id}",
            sep="\t"
        )

        print(
            seq_name,
            source_name,
            "mRNA",
            exon_pos_list[0][0] + 1,
            exon_pos_list[-1][1],
            ".",
            strand,
            ".",
            f"ID={trans_id};Parent={gene_id}",
            sep="\t",
        )

        five_utr_pos_list = left_pos_list
        three_utr_pos_list = right_pos_list
        if strand == "-":
            exon_pos_list.sort(key=lambda x: x[0], reverse=True)
            center_pos_list.sort(key=lambda x: x[0], reverse=True)
            left_pos_list.sort(key=lambda x: x[0], reverse=True)
            right_pos_list.sort(key=lambda x: x[0], reverse=True)
            five_utr_pos_list = right_pos_list
            three_utr_pos_list = left_pos_list

        for i in range(0, len(exon_pos_list)):
            print(
                seq_name,
                source_name,
                "exon",
                exon_pos_list[i][0] + 1,
                exon_pos_list[i][1],
                ".",
                strand,
                ".",
                f"ID={trans_id}.exon{i + 1};Parent={trans_id}",
                sep="\t",
            )

        for i in range(0, len(five_utr_pos_list)):
            print(
                seq_name,
                source_name,
                "five_prime_UTR",
                five_utr_pos_list[i][0] + 1,
                five_utr_pos_list[i][1],
                ".",
                strand,
                ".",
                f"ID={trans_id}.5utr{i + 1};Parent={trans_id}",
                sep="\t",
            )

        for i in range(0, len(center_pos_list)):
            print(
                seq_name,
                source_name,
                "CDS",
                center_pos_list[i][0] + 1,
                center_pos_list[i][1],
                ".",
                strand,
                ".",
                f"ID={trans_id}.cds{i + 1};Parent={trans_id}",
                sep="\t",
            )

        for i in range(0, len(three_utr_pos_list)):
            print(
                seq_name,
                source_name,
                "three_prime_UTR",
                three_utr_pos_list[i][0] + 1,
                three_utr_pos_list[i][1],
                ".",
                strand,
                ".",
                f"ID={trans_id}.3utr{i + 1};Parent={trans_id}",
                sep="\t",
            )
