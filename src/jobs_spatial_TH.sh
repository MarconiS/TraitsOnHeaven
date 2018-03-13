for datafile in "$@"
do
  filename="${datafile##*/}"
  echo $filename
  sbatch siteScale.SLURM $filename OSBS
done
