#!/bin/bash
#interpreter used to execute the script

#“#SBATCH” directives that convey submission options:
files=("s1" "s2")
j=0
for i in $(find /scratch/gthurber_root/gthurber0/marcase/NovaSeq5May22/fastqs_5829-MC  -name '*.fastq'); do
	j=$((j+1))
	job_dir="${i}"
	echo "job_dir: ${job_dir}"
	job_file="$(basename $job_dir)"
	echo "job_file: ${job_file}"

	echo "#!/bin/bash


#SBATCH --job-name=${job_file}
#SBATCH --mail-user=marcase@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64000m
#SBATCH --time=8:00:00
#SBATCH --account=gthurber0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j-${job_file//.fastq/}_merge.log

# The application(s) to execute along with its input arguments and options:
module load python3.8-anaconda/2021.05
python NGS_pattern_finder_combined.py ${job_dir}" > ${job_dir//.fastq}.sh
	sbatch ${job_dir//.fastq/}.sh
done
