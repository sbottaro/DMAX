import numpy as np
import mdtraj as md
import sys
import itertools as its

traj = sys.argv[1]
if(len(sys.argv)==3):
    top = sys.argv[2]
    tt = md.load(traj,top=top)
else:
    tt = md.load(traj)

#md.compute_rdf()    
#pairs = [[[i,j] for j in range(i+1,len(tt.atoms))]  ]
ll= tt.topology.n_atoms
pairs = list(its.combinations(range(ll),2))
#print(len(pairs))
dists = md.compute_distances(tt,pairs)
#print(dists.shape)
imax = np.argmax(dists[0,:])
atoms = [str(at) for at in tt.topology.atoms]
#print(dists[0,imax],pairs[imax],atoms[pairs[imax][0]],atoms[pairs[imax][1]])
print("%s DMAX=%-8.3f[nm] (%s/%s)" % (traj,dists[0,imax],atoms[pairs[imax][0]],atoms[pairs[imax][1]]))
