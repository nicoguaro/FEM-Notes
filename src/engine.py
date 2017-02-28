"""
PROGRAM SOLIDS
--------------

Computes the displacement solution for an assembly
of 1D spring elements under point loads.

Created by Juan Gomez as part of the courses:

IC0283 COMPUTATIONAL MODELLING
IC0602 INTRODUCTION TO THE FINITE ELEMENT METHOD
Universidad EAFIT
Departamento de Ingenieria Civil
Last updated February 2017
"""
import numpy as np
import starter as sta
import assembler as ass
from datetime import datetime
start_time = datetime.now()
"""
   PRE-PROCESSING
"""
nodes , mats , elements , loads = sta.readin()
ne , nn , nm , nl = sta.proini(nodes , mats , elements , loads)
DME , IBC , neq = ass.DME(nn , ne , nodes , elements)
"""
   ASSEMBLY
"""
KG = np.zeros([neq, neq])
for i in range(ne):
    kloc , ndof  = ass.retriever(elements , mats  , nodes , i)
    KG = ass.assembler(KG , neq , kloc , ndof , DME , i)
RHSG = ass.loadasem(loads, IBC, neq, nl)
"""
   SYSTEM SOLUTION
"""
UG = np.linalg.solve(KG, RHSG)
if not(np.allclose(np.dot(KG, UG), RHSG)):
    print("The system is not in equilibrium!")
end_time = datetime.now()
print('Duration for system solution: {}'.format(end_time - start_time))



























