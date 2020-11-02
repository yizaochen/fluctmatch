#!/bin/bash
# Program:
#       This program sets up dcd and pdb of central base-pairs

# Ad hocs
rootfolder=$1
host=$2
type_na=$3

if [ $type_na == "arna+arna" ]; then
    na=arna
elif [ $type_na == "bdna+bdna" ]; then
    na=bdna
fi

gmx=/home/yizaochen/opt/gromacs/bin/gmx
inputfolder=${rootfolder}/${host}/${type_na}/input/allatoms
inp_xtc=${inputfolder}/${type_na}.all.xtc
inp_gro=${inputfolder}/${type_na}.npt4.all.gro

# gro to pdb
out_pdb=${inputfolder}/${type_na}.npt4.all.pdb
${gmx} editconf -f ${inp_gro} -o ${out_pdb}

# make index
ndx=${inputfolder}/${type_na}.ndx
${gmx} make_ndx -f ${inp_gro} -o ${ndx}

# make central pdb
central_pdb=${inputfolder}/${type_na}.central.pdb
${gmx} editconf -f ${inp_gro} -o ${central_pdb} -n ${ndx}

# make central xtc
central_xtc=${inputfolder}/${type_na}.central.xtc
${gmx} trjconv -s ${inp_gro} -f ${inp_xtc} -o ${central_xtc} -n ${ndx}

# cp central pdb to pdb1 pdb2
pdb1=${inputfolder}/${na}1.central.pdb
pdb2=${inputfolder}/${na}2.central.pdb
cp $central_pdb $pdb1
cp $central_pdb $pdb2
