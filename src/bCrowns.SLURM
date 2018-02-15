#!/bin/bash
#SBATCH --job-name=getCrData    # Job name
#SBATCH --mail-type=FAIL,END               # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite -b 

#SBATCH --ntasks=1                 # Number of MPI ranks
#SBATCH --cpus-per-task=15            # Number of cores per MPI rank 
#SBATCH --mem=160000
#SBATCH --time=90:00:00       #Time limit hrs:min:sec
#SBATCH --output=getCr.out   # Standard output and error log
#SBATCH --error getCr.err


pwd; hostname; date
ml gcc/6.3.0 R/3.4.0 python gdal/2.2.1 
echo "getCrData"
cd /ufrc/ewhite/s.marconi/Marconi2018
Rscript --default-packages=stats,utils,datasets,base,graphics,grDevices,methods --vanilla spatialCore.R $1
echo "done"

date
