grep '.' results.txt | awk '$NF < 1.5' | awk '{print $NF}' > rmsd.list
grep '.' results.txt | awk '$NF < 1.5' | awk -F '.' '{print $1}' > index.list
#for i in `grep '.' results.txt | awk '$NF < 1.5' | awk -F '.' '{print $1}'` ; do grep SMILES poses.mol2 | awk '{print $NF}' | head -${i} | tail -1 ; done > smiles.list
for i in `grep '.' results.txt | awk '$NF < 1.5' | awk -F '.' '{print $1}'` ; do grep Name poses.mol2 | grep -v 'Long' | awk '{print $NF}' | head -${i} | tail -1 | awk -F '.' '{print $1}' ; done > names.list
for name in `cat names.list ` ; do pwd | awk -F 'CYS' '{print $NF}' | awk -v "var=${name}" -F '/' '{print "cat /home/danielza/Work/Current/Covalentizer/databases/PDB_Final/Ligands/Prepared/" $2 "/linkers_index.smi | awk '\''$2==\""var"\"'\'' | awk '\''{print $1}'\''"}' ; done > tmp.sh ; chmod 777 tmp.sh ; ./tmp.sh > smiles.list
rm tmp.sh
paste smiles.list names.list rmsd.list index.list | sort -k3,3 -n | awk '!x[$1]++' | cat -n > data.txt
python $COVALENTIZER/Scripts/ExtractPoses.py poses.mol2 index.list
for i in `wc -l data.txt | awk '$1 > 0' | awk '{print $1}'` ; do mkdir web_files/ ; done
cat data.txt | awk '{print "cp poses/"$5".mol2 web_files/"$1".mol2 ; echo \""$2"\" > web_files/"$1".smi ; echo "$4" > web_files/"$1".rmsd"}' > make_web_folder.sh ; chmod 777 make_web_folder.sh ; ./make_web_folder.sh
if [ -d "web_files/" ]; then
    cp ../rec.pdb web_files/
    cp ../xtal-lig.pdb web_files/
    head -1 ../../../res.txt | awk '{print $2}' > web_files/cys_position.txt
    head -1 ../../../res.txt | awk '{print $3}' > web_files/chain.txt
    cd web_files/
    for i in *.smi ; do python $COVALENTIZER/Scripts/reform_elec.py $COVALENTIZER/Scripts/reform_elec.re $i tmp.smi ; cat tmp.smi | awk '{print $1}' > $i ; python $COVALENTIZER/Scripts/draw_smile.py $i ; done ; rm tmp.smi
    cd ../
fi
