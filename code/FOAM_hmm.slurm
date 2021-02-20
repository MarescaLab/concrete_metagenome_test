#!/bin/bash
#SBATCH --job-name=FOAM_hmm
#SBATCH --chdir=/work/akiledal/concrete_metagenome_test/data/FOAM
#SBATCH --mem=500G
#SBATCH -c 64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=akiledal@udel.edu
#SBATCH --time=20-0

source ~/.bashrc

date
hostname
start=`date +%s`

###Script starts here###

conda activate hmmer

which hmmer

hmmsearch -h

hmmsearch \
	--cpu 64 \
	--tblout s3_sixframe.txt \
	--domT 14 \
	-o s3_sixframe.hm \
	../reference/FOAM-hmm_rel1.hmm \
	six_frame_S3_f.fasta

tail -n+4 s3_sixframe.txt | sed 's/ * / /g' | cut -f 1-10 -d " " > s3_sixframe_mod.txt

###Script ends here###
echo "Duration: $((($(date +%s)-$start)/60)) minutes"