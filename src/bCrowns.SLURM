#!/bin/bash
#SBATCH --job-name=getCrData    # Job name
#SBATCH --mail-type=END               # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --ntasks=1                 # Number of MPI ranks
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank
#SBATCH --mem=4000
#SBATCH --time=00:59:00       #Time limit hrs:min:sec
#SBATCH --output=itc.out   # Standard output and error log
#SBATCH --error=itc.err


pwd; hostname; date
ml gcc/6.3.0 R/3.4.0 python gdal/2.2.1
echo "spatial_core_2017"
echo $1

Rscript --default-packages=stats,utils,datasets,base,graphics,grDevices,methods --vanilla ./src/getData.R $1 $2
echo "done"

date
