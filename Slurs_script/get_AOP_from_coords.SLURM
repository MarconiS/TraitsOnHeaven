#!/bin/bash
#SBATCH --job-name=extract_traits    # Job name
#SBATCH --mail-type=BEGIN,FAIL,END               # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --ntasks=1                 # Number of MPI ranks
#SBATCH --cpus-per-task=32            # Number of cores per MPI rank
#SBATCH --mem=250000

#SBATCH --time=36:00:00       #Time limit hrs:min:sec
#SBATCH --output=extract_traits.out   # Standard output and error log
#SBATCH --error=extract_traits.err


pwd; hostname; date
ml R/3.4.3 conda
source activate quetzal3
echo "grepping data from tiles"
Rscript --default-packages=stats,utils,datasets,base,methods  ../AOP_from_coords/src/main.R $1
echo "done"

date

