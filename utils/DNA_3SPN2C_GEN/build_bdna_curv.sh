#!/bin/sh

# Calls to


if [ $# -ne 1 ]; then
    echo "Usage: $0 <sequence file>"
    exit 1
fi

ICNF=./ICNF/icnf.exe
UTILS=.

echo "Making parameter file"
./make_bp_params.py $1

echo "Running X3DNA"
export PATH=$PATH:~/Workspace/x3dna-v2.1/bin
export X3DNA=~/Workspace/x3dna-v2.1
x3dna_utils cp_std BDNA
rebuild -atomic bp_step.par atomistic.pdb
mkdir basis
mv Atomic* basis/

echo "Mapping to CG coordinates"
./pdb2cg_dna.py atomistic.pdb
mv in00_conf.xyz bdna_curv.xyz

echo "Building LAMMPS input file"
${ICNF} $1 1 1 . 0

echo "Replacing atoms in configuration file"
./replace_atoms.sh conf_lammps.in dna_conf.in bdna_curv_conf.in

echo "Making list files"
./make_list_files.py bdna_curv_conf.in

echo "Making ninfo files"
./make_ninfo_files.py $1
