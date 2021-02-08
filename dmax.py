import numpy as np
import mdtraj as md
import sys
import itertools as its
from scipy.spatial.distance import pdist,squareform
import argparse

def calc(name,tt):

    ll= tt.topology.n_atoms
    atoms = [str(at) for at in tt.topology.atoms]
    
    for k in range(len(tt)):
        # first calculate distance between center of mass
        coms = []
        ats_idxs = []
        for res in tt.topology.residues:
            ii = [at.index for at in res.atoms]
            coms.append(np.average(tt.xyz[k,ii],axis=0))
            ats_idxs.append(ii)
        if(len(coms)==1):
            print("ERROR: FILE %s contains only one residue" % name)
            break
        com_dist = pdist(coms)
	
        # calculate pairwise distances only between top 2% 
        cut = int(len(com_dist)*0.02)
        com_dist_max = np.argsort(com_dist)[-cut:]

        triu_idx = np.triu_indices(len(coms),1)
        far_residues = [[triu_idx[0][k],triu_idx[1][k]] for k in com_dist_max]
        pairs = []
        labels = []

        for j in far_residues:
            pairs.extend([(x,y) for x in ats_idxs[j[0]] for y in ats_idxs[j[1]]])
            labels.extend(["%s/%s" % (atoms[x],atoms[y]) for x in ats_idxs[j[0]] for y in ats_idxs[j[1]]] )
	
            
            #print("# N residuese=%d; residue pairs=%d; Pruned=%d; Total calc=%d/%d" % (len(coms),len(com_dist),cut,len(pairs)+len(com_dist),ll*(ll-1)/2))
            # calculate distances between atoms in distant residues
        dists = md.compute_distances(tt[k],pairs)
        imax = np.argmax(dists[0,:])
        print("%10s %8d %-8.4f %s" % (name,k,dists[0,imax],labels[imax]))
	

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser("Calculate Dmax of PDBs or trajectory")
    parser.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser.add_argument("--top", dest="top",help="Topology file",required=False)
    args = parser.parse_args()
    #print(args.pdbs)
    if(args.pdbs!=None):
        assert(args.trj==None)
        assert(args.top==None)
        for pdb in args.pdbs:
            calc(pdb.split(".pdb")[0],md.load(pdb))
    else:
        assert(args.pdbs==None)
        calc(md.load(args.trj.split(".")[0],args.trj,top=args.top))
