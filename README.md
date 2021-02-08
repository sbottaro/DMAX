# DMAX
Calculate the maximum distance between atoms in a PDB file or trajectory.

- Mdtraj, scipy and numpy are required. 
- USAGE: python dmax.py --pdb file1.pdb file2.pdb ... fileN.pdb 
- USAGE: python dmax.py --top topology_file --traj trajectory_file 

dmax_full.py is very slow and should not be used unless dmax.py gives strange results. 

