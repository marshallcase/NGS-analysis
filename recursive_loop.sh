#!/bin/bash
#interpreter used to execute the script

#“#SBATCH” directives that convey submission options:
files=("s1" "s2")
j=0
for i in $(find /scratch/gthurber_root/gthurber0/marcase/NovaSeq28Feb22/fastqs_5243-MC  -name '*_R1_*.fastq'); do
	j=$((j+1))
	job_dir="${i//.fastq/}"
	echo "${job_dir}.sh"
	job_file="$(basename $job_dir)"
	#echo "${job_dir}.fastq"

	echo "#!/bin/bash


#SBATCH --job-name=${job_file}
#SBATCH --mail-user=marcase@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64000m 
#SBATCH --time=4:00:00
#SBATCH --account=gthurber0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j-${job_file}_n1.log

# The application(s) to execute along with its input arguments and options:
module load python3.8-anaconda/2021.05 
python /scratch/gthurber_root/gthurber0/marcase/NovaSeq28Feb22/NGS_SeqIO_lakes.py ${job_dir}.fastq ${job_dir//_R1_/_R2_}.fastq" > ${job_dir}_n1.sh
	sbatch ${job_dir}_n1.sh

done
