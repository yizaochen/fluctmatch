#!/bin/bash
# Program:
#       This program convert gromac's output pdb into pdb which charmm can process 

rootfolder=$1
host=$2
type_na=$3
strand_num=$4

if [ $type_na == "arna+arna" ]; then
        na=arna
elif [ $type_na == "bdna+bdna" ]; then
        na=bdna
fi

inputfolder=${rootfolder}/${host}/${type_na}/input/allatoms
outputfolder=${rootfolder}/${host}/${type_na}/make_crd
head=/home/yizaochen/codes/na_mechanics/shell_scripts/head.pdb

SYSTEM=${na}${strand_num}
cp ${head} ${outputfolder}/${SYSTEM}.1.pdb

#grep "ATOM   " $fid0.pdb | grep " A " | sed  's/           [A-Z]/      '$fid'/' >> $fid.0.pdb
grep "ATOM  " ${inputfolder}/${SYSTEM}.central.pdb |\
sed 's/RU5P/URA /g'|\
sed 's/DT5P/THY /g'|\
sed 's/DT5/THY/g'|\
sed 's/DT3/THY/g'|\
sed 's/DA5/ADE/g'|\
sed 's/DA3/ADE/g'|\
sed 's/DG5/GUA/g'|\
sed 's/DG3/GUA/g'|\
sed 's/DC5/CYT/g'|\
sed 's/DC3/CYT/g'|\
sed 's/RA5/ADE/g'|\
sed 's/RA3/ADE/g'|\
sed 's/RU5/URA/g'|\
sed 's/RU3/URA/g'|\
sed 's/RC5/CYT/g'|\
sed 's/RC3/CYT/g'|\
sed 's/RG5/GUA/g'|\
sed 's/RG3/GUA/g'|\
sed 's/ DT/THY/g'|\
sed 's/ DG/GUA/g'|\
sed 's/ DC/CYT/g'|\
sed 's/ DA/ADE/g'|\
sed 's/ RA/ADE/g'|\
sed 's/ RU/URA/g'|\
sed 's/ RG/GUA/g'|\
sed 's/ RC/CYT/g'|\
sed 's/1H5'"'"'/ H5'"'"'/g'|\
sed 's/2H5'"'"'/H5'"'"''"'"'/g'|\
sed 's/1H2'"'"'/ H2'"'"'/g'|\
sed 's/2H2'"'"'/H2'"'"''"'"'/g'|\
sed 's/2HO'"'"'/H2'"'"''"'"'/g'|\
sed 's/1H6/H61/g' | \
sed 's/2H6/H62/g' | \
sed 's/1H4/H41/g' | \
sed 's/2H4/H42/g' \
>> ${outputfolder}/$SYSTEM.1.pdb
echo 'END' >> ${outputfolder}/$SYSTEM.1.pdb
