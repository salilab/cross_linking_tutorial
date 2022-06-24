#!/usr/bin/env python3


from __future__ import print_function

import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.io
import IMP.pmi.io.crosslink

xldb='''Protein 1,Protein 2,Residue 1,Residue 2,UniqueID,Score
ProtA,ProtB,1,10,1,1.0
ProtA,ProtB,1,11,1,2.0
ProtA,ProtB,1,21,2,2.0
'''

with open('xlinks.csv', 'w') as xlf:
    xlf.write(xldb)

cldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_residue2_key("Residue 2")
cldbkc.set_unique_id_key("UniqueID")
cldbkc.set_id_score_key("Score")

cldb = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("xlinks.csv")

print(cldb)

xldb='''Protein 1,Protein 2,Residue 1,Residue 2,UniqueID,Score
ProtA,ProtB,1,10,1,1.0
ProtA,ProtB,1,11,1,2.0
ProtB,ProtA,21,1,2,2.0
ProtA,ProtA,1,2,3,3.0
'''

with open('xlinks.csv', 'w') as xlf:
    xlf.write(xldb)

cldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_residue2_key("Residue 2")
cldbkc.set_unique_id_key("UniqueID")
cldbkc.set_id_score_key("Score")

cldb = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("xlinks.csv")

from IMP.pmi.io.crosslink import FilterOperator as FO
import operator

fo1 = FO(cldb.protein1_key, operator.eq, "ProtA")
cldb.set_value(cldb.protein1_key, "ProtA.1", fo1)
fo2 = FO(cldb.protein2_key, operator.eq, "ProtA")
cldb.set_value(cldb.protein2_key, "ProtA.1", fo2)

cldb.clone_protein("ProtA.1", "ProtA.2")

print(cldb)
