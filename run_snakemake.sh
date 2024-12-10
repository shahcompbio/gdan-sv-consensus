cluster_fmt="sbatch --partition=componc_cpu --cpus-per-task={threads} --mem={resources.mem_mb} --job-name={rule}.{wildcards} --error=logs/{rule}/{rule}.{wildcards}.%j.err --output=logs/{rule}/{rule}.{wildcards}.out --time=24:00:00"
cmd="snakemake --executor cluster-generic"
cmd="$cmd --configfile config.yaml"
cmd="$cmd --cluster-generic-submit-cmd \"$cluster_fmt\""
cmd="$cmd --profile profile/"
cmd="$cmd --singularity-args \"--bind /data1 --bind /home\""
cmd="$cmd --dryrun"
cmd="$cmd --rulegraph"
cmd="$cmd -p"
echo $cmd
eval $cmd
