#!/usr/bin/env python


class GFF3Tree:
    def __init__(self):
        self.parent_dict = dict()

    def set_relations(self, fn):
        with open(fn) as fin:
            for ln in fin:
                if len(ln) == 0 or ln[0] == "#":
                    continue
                f = ln.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attr_str = f[8]

                feat_id = ""
                parent_id = ""
                for attr in attr_str.split(";"):
                    attr_pair = attr.split("=")
                    if len(attr_pair) != 2:
                        continue
                    if attr_pair[0] == "ID":
                        feat_id = attr_pair[1]
                    elif attr_pair[0] == "Parent":
                        parent_id = attr_pair[1]
                if feat_id != "" and parent_id != "":
                    self.parent_dict[feat_id] = parent_id

    def get_parent_id(self, feat_id):
        if feat_id in self.parent_dict:
            return self.parent_dict[feat_id]
        else:
            return ""

    def get_root_id(self, feat_id):
        while feat_id in self.parent_dict:
            feat_id = self.parent_dict[feat_id]
        return feat_id


def gff3_attr_get_id(attr_str):
    feat_id = ""
    for attr in attr_str.split(";"):
        attr_pair = attr.split("=")
        if len(attr_pair) == 2 and attr_pair[0] == "ID":
            feat_id = attr_pair[1]
            break
    return feat_id
