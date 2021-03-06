# DETAILS FROM THE CONFIG FILE 
run_he_script="/full/path/to/experiment_name/H276_GABA_runHyper.sh" # if needed insert the proper path before the script name
bwa_aln_soft="full/path/to/bwa/0.7.15/bwa-0.7.15/bwa" # if needed insert the proper path before the tool name
bwa_mem_soft="full/path/to/bwa/0.7.15/bwa-0.7.15/bwa" # if needed insert the proper path before the tool name
SamToFastq_soft="full/path/to/samtools/1.9/bin/samtools" # if needed insert the proper path before the tool name
SamTools_soft="full/path/to/samtools/1.9/bin/samtools" # if needed insert the proper path before the tool name

# MODIFY PER SAMPLE
echo -e "/full/path/to/H276_GABA.Unmapped.out.mate1.gz\tH276_GABA_out\n/full/path/to/H276_GABA.Unmapped.out.mate2.gz\tH276_GABA_out" > /full/path/to/experiment_name/H276_GABA_file_list
source_file="/full/path/to/experiment_name/H276_GABA_file_list" # path+full_name of the input files to run: fastq_file_path+name	/TAB/ out_prefix (if the input files are of paired-end reads, each of the paired files should appear in separate line)
dir_pre="/full/path/to/experiment_name/H276_GABA" # path+full_name of the output directory

# HARD CODED
genome_bwa_ind="/full/path/to/experiment_name/HE_GenomeIndex/GRCh38.primary_assembly.genome" # path+prefix of the index genome expected 5 files like: prefix.amb, prefix.ann, prefix.bwt, prefix.pac, prefix.sa
genome_trans_bwa_ind="/full/path/to/experiment_name/HE_GenomeIndex_Trans/GRCh38.primary_assembly.genome" # path+prefix of the index transformed genome: for each of the 6 transformed (a2c a2g a2t g2c t2c t2g)=>tt: 5 files: prefix.tt.amb, prefix.tt.ann, prefix.tt.bwt, prefix.tt.pac, prefix.tt.sa + 1 fasta file prefix.tt.fa => tot 36 files
genome_fasta="/full/path/to/experiment_name/HE_GenomeIndex/GRCh38.primary_assembly.genome.fa" # path+full_name of the fasta file of the original genome
Q="33" # PHRED score offset of base quality (33 or 64)
PE="1" # 0 => single end, 1 => paired end
GAP=500000 # gap max size between the pairs
bwa_run="1" # 0 => pipeline will not run BWA mapping if analyseMM file exists, 1 => first run BWA mapping of the source files  
ue_detect_args="0.05 0.6 30 0.6 0.1 0.8 0.2" # args meaning: -Min of edit sites at Ultra-Edit read -Min fraction of edit sites/mm sites -Min sequence quality for counting editing event -Max fraction of same letter in cluster -Min of cluster length -Max initiate index of cluster -Min ending index of cluster

################################################################################################

sh $run_he_script $genome_bwa_ind $genome_trans_bwa_ind $genome_fasta $Q $PE $GAP $dir_pre $bwa_run $ue_detect_args $source_file $bwa_aln_soft $bwa_mem_soft $SamToFastq_soft $SamTools_soft
