#!/usr/bin/env python3


from __future__ import print_function

import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
import ihm.cross_linkers

m = IMP.Model()

s = IMP.pmi.topology.System(m)
 
st1 = s.create_state()
protA = st1.create_molecule("ProtA", "G" * 10, "A")
protA.add_representation(resolutions=[10])
protB = st1.create_molecule("ProtB", "G" * 30, "B")
protB.add_representation(resolutions=[10])
hier = s.build()

beads = IMP.atom.Selection(hier).get_selected_particles()
print(beads)
xyzs = [IMP.core.XYZ(b) for b in beads if IMP.core.XYZ.get_is_setup(b)]
xyzs[0].set_coordinates((0,0,0))
xyzs[1].set_coordinates((-40,0,0))
xyzs[2].set_coordinates((0,0,0))
xyzs[3].set_coordinates((40,0,0))

xldb='''Protein 1,Protein 2,Residue 1,Residue 2,UniqueID,Score
ProtA,ProtB,1,10,1,1.0
ProtA,ProtB,1,11,2,2.0
ProtA,ProtB,1,21,3,2.0
'''
with open('xlinks.csv', 'w') as fh:
    fh.write(xldb)

cldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_residue2_key("Residue 2")

# the unique_id_key and id_score_key are optional,
# and they add features that will be explained below

cldbkc.set_unique_id_key("UniqueID")
cldbkc.set_id_score_key("Score")

cldb = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("xlinks.csv")

print(cldb)

xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=hier, database=cldb, length=21.0, slope=0.0,
    resolution=1.0, label="XL", linker=ihm.cross_linkers.dss)

print(xl.rs.unprotected_evaluate(None))

sel = IMP.atom.Selection(hier, molecule="ProtA")
pA, = sel.get_selected_particles()

scores = []
xs = []
for i in range(-100, 100):
    xs.append(float(i))
    IMP.core.XYZ(pA).set_coordinates((i, 0, 0))
    scores.append(xl.rs.unprotected_evaluate(None))

import pylab

pylab.plot(xs, scores)

from IMP.pmi.io.crosslink import FilterOperator
import operator

fo = FilterOperator(cldb.unique_id_key, operator.eq, "2")
fcldb = cldb.filter(fo)

print(fcldb)

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=hier, database=fcldb, length=21.0, slope=0.0,
    resolution=1.0, label="XL", linker=ihm.cross_linkers.dss)

scores=[]
xs=[]
for i in range(-100, 100):
    xs.append(float(i))
    IMP.core.XYZ(pA).set_coordinates((i, 0, 0))
    scores.append(xl1.rs.unprotected_evaluate(None))

pylab.plot(xs, scores)

sigma = xl1.sigma_dictionary["SIGMA"][0]

scores_list = []
xs_list = []
for s in range(1, 20):
    scores = []
    xs = []
    sigma.set_scale(float(s))
    for i in range(-100, 100):
        xs.append(float(i))
        IMP.core.XYZ(pA).set_coordinates((i, 0, 0))
        scores.append(xl1.rs.unprotected_evaluate(None))
    scores_list.append(scores)
    xs_list.append(xs)

for xs, scores in zip(xs_list, scores_list):
    pylab.plot(xs, scores)

sigma.set_scale(11)
psi = xl1.psi_dictionary["PSI"][0]

import numpy as np

scores_list = []
xs_list = []
for s in np.linspace(0.01, 0.5, 10):
    scores = []
    xs = []
    psi.set_scale(float(s))
    for i in range(-100, 100):
        xs.append(float(i))
        IMP.core.XYZ(pA).set_coordinates((i, 0, 0))
        scores.append(xl1.rs.unprotected_evaluate(None))
    scores_list.append(scores)
    xs_list.append(xs)

for xs, scores in zip(xs_list, scores_list):
    pylab.plot(xs, scores)

sigma = xl.sigma_dictionary["SIGMA"][0]

scores_list = []
xs_list = []
for s in range(1, 20):
    scores = []
    xs = []
    sigma.set_scale(float(s))
    for i in range(-100, 100):
        xs.append(float(i))
        IMP.core.XYZ(pA).set_coordinates((i, 0, 0))
        scores.append(xl.rs.unprotected_evaluate(None))
    scores_list.append(scores)
    xs_list.append(xs)

for xs, scores in zip(xs_list, scores_list):
    pylab.plot(xs, scores)

sigma.set_scale(11)
psi = xl.psi_dictionary["PSI"][0]

import numpy as np

scores_list = []
xs_list = []
for s in np.linspace(0.01, 0.5, 10):
    scores = []
    xs = []
    psi.set_scale(float(s))
    for i in range(-100, 100):
        xs.append(float(i))
        IMP.core.XYZ(pA).set_coordinates((i, 0, 0))
        scores.append(xl.rs.unprotected_evaluate(None))
    scores_list.append(scores)
    xs_list.append(xs)

for xs, scores in zip(xs_list, scores_list):
    pylab.plot(xs, scores)

import math
import numpy as np

sel = IMP.atom.Selection(hier, molecule="ProtB")
pB1, pB2, pB3 = sel.get_selected_particles()

IMP.core.XYZ(pA).set_coordinates((0,0,0))
IMP.core.XYZ(pB1).set_coordinates((0,20,0))
IMP.core.XYZ(pB2).set_coordinates((-20*math.sqrt(3)/2,-20/2,0))
IMP.core.XYZ(pB3).set_coordinates((20*math.sqrt(3)/2,-20/2,0))

scores = []
psis = []
sigmas = []
for p in np.linspace(0.01, 0.5, 100):
    psi.set_scale(p)
    for s in np.linspace(1, 40, 50):
        psis.append(p)
        sigmas.append(s)
        sigma.set_scale(s)
        scores.append(xl.rs.unprotected_evaluate(None))

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.scatter(psis, sigmas, c=scores, s=30, edgecolor=[])
plt.show()

IMP.core.XYZ(pA).set_coordinates((100, 100, 100))
scores = []
psis = []
sigmas = []
for p in np.linspace(0.01, 0.5, 100):
    psi.set_scale(p)
    for s in np.linspace(1, 40, 50):
        psis.append(p)
        sigmas.append(s)
        sigma.set_scale(s)
        scores.append(xl.rs.unprotected_evaluate(None))

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.scatter(psis, sigmas, c=scores, s=30, edgecolor=[])
plt.show()
