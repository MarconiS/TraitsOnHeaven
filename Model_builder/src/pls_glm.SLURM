#!/bin/bash
#SBATCH --job-name=traitsBuild    # Job name
#SBATCH --mail-type=BEGIN,FAIL,END               # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --ntasks=1                 # Number of MPI ranks
#SBATCH --cpus-per-task=64            # Number of cores per MPI rank
#SBATCH --mem=24000

#SBATCH --time=24:05:00       #Time limit hrs:min:sec
#SBATCH --output=plamBuild.out   # Standard output and error log
#SBATCH --error=plmaBuild.err


pwd; hostname; date
ml gcc/6.3.0 R/3.4.0 python gdal/2.2.1
echo "Building model"
Rscript --default-packages=stats,utils,datasets,base,methods --vanilla ./Model_builder/src/pCore.R $1 $2
echo "done"

date
