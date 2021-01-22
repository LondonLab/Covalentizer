cat extr*sort.txt | awk '$NF < 10' > new_extract.txt ; /work/londonlab/git_dock/DOCK/analysis/getposes.py -f new_extract.txt -x 0 ; python /home/danielza/CovaLib/Scripts/SeperatePoses.py poses.mol2
