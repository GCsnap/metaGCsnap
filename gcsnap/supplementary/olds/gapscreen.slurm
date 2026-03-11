#!/bin/bash                  
#The previous line is mandatory

#SBATCH --job-name=gap_screen     #Name of your job
#SBATCH --cpus-per-task=60    #Number of cores to reserve
#SBATCH --mem=150G     #Amount of RAM/core to reserve
#SBATCH --time=24:00:00      #Maximum allocated time
#SBATCH --qos=1week        #Selected queue to allocate your job
#SBATCH --output=log/gap_screen.o%j   #Path and name to the file for the STDOUT
#SBATCH --error=log/gap_screen.e%j    #Path and name to the file for the STDERR

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ribosome

python gapscreen.py --threads 6 --batch-size 10