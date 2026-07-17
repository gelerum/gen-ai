#!/usr/bin/env python3
import numpy as np
from goatools.obo_parser import GODag

GO_DAG = GODag("../go-basic.obo")

TRANSPORT_ROOTS = {
    "GO:0006810",  # transport
    "GO:0005215",  # transporter activity
    "GO:0055085",  # transmembrane transport
    "GO:0015031",  # protein transport
    "GO:0008565",  # protein transporter activity
}

transport_terms = set()

for go_id in TRANSPORT_ROOTS:

    if go_id not in GO_DAG:
        continue

    term = GO_DAG[go_id]

    transport_terms.add(go_id)

    for child in term.get_all_children():
        transport_terms.add(child)


def is_transport_protein(go_annotations):
    if type(go_annotations) != np.ndarray:
        return None

    for go in go_annotations:

        go_id = go["id"]

        if go_id in transport_terms:
            return True

    return False