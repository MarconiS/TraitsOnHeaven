for datafile in "$@"
do
  filename="${datafile##*/}"
  echo $filename
  sbatch bCrowns.SLURM $filename
done