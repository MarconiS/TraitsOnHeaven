#!/bin/bash
#SBATCH --job-name=traitsBuild    # Job name
#SBATCH --mail-type=BEGIN,FAIL,END               # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --ntasks=1                 # Number of MPI ranks
#SBATCH --cpus-per-task=64            # Number of cores per MPI rank
#SBATCH --mem=240000

#SBATCH --time=36:05:00       #Time limit hrs:min:sec
#SBATCH --output=plamBuild.out   # Standard output and error log
#SBATCH --error=plmaBuild.err


pwd; hostname; date
ml R/3.4.3 conda
source activate quetzal3

echo "Building model"
Rscript ./AOP_from_coords/src/main.R
echo "done"

date
