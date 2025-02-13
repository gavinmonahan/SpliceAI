#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --job-name=spliceai-gpu
#SBATCH --partition=gpu-dev
#SBATCH --account=pawsey0933-gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --export=NONE
#SBATCH --mail-user=gavin.monahan@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out

module unload gcc/12.2.0
module swap pawseyenv/2024.05 pawseyenv/2023.08
module load gcc/12.2.0
module load tensorflow/rocm5.6-tf2.12

bash

cd $MYSOFTWARE/manual/software/pythonEnvironments/tensorflowContainer-environments

source spliceai2/bin/activate

/home/gmonahan/.local/bin/spliceai -R /scratch/references/human/broad_hg19/Homo_sapiens_assembly19.fasta -A grch37 -I $MYSOFTWARE/T2T/SpliceAI/examples/input.vcf -O $MYSOFTWARE/T2T/SpliceAI/examples/test_out.vcf
