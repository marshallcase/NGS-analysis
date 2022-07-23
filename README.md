# NGS-analysis
analyze bacterial surface display NGS data sets 

step 1)

Use NGMerge (https://github.com/jsh58/NGmerge) to align and merge forward and reverse .fastq files

One at a time:

/NGmerge  -1 sample_R1.fastq.gz  -2 sample_R2.fastq.gz  -o sample_merged.fastq.gz

Or use recursive_loop.gz.sh to automate it

step 2)

identify correctly expressed and read bacterial clones

One at a time:

NGS_pattern_finder_combined.py

Or use recursive_loop_aligned.sh to automate it

step 3)

combine individual read files into a single dataframe object for analysis

Use NGS_merge.py and merge.sh

step 4a)

read counts per round of sorting: analyze data using analysis_script.py

step 4b)

analyze SORTCERY data using reads_plots.py
analyze gate scores using gate_score_analysis.py

step 5)

machine learning - explore models and performance using ml_script.py

step 6)

perform ILP - ILP_order2.py


