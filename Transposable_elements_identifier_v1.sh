#!/usr/local_rwth/bin/bash

#SBATCH --job-name=Transposable_elements_identifier_jiangzhao
#SBATCH --output=Transposable_elements_identifier_jiangzhao.txt
#SBATCH --error=Transposable_elements_identifier_jiangzhao.error
#SBATCH --nodes=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=12:00:00
#SBATCH -A rwth0146

#SBATCH --mail-type=END
#SBATCH --mail-user=jqian@bio1.rwth-aachen.de

export PATH=$PATH:/home/ur376715/miniconda3/bin


### EDTA can be installed by conda https://github.com/oushujun/EDTA
### RepeatMasker was included inside EDTA package

### the directory of repeatmasker library. Later the species specific Repbase (2018_10) will be integrated here

### EDTA results will be stored here
work_dir=/work/ur376715/Second_time_EDTA_result

RM_dir=/home/ur376715/miniconda3/envs/EDTA/share/RepeatMasker/Libraries

Genome_dir=/home/ur376715/Bgh/Bgh_seq/Genome/bgh_dh14_v4.fa

Genome_CDS_dir=/home/ur376715/Bgh/Bgh_seq/Genome/bgh_dh14_v4-2.cds.fa

Genome_BED_dir=/home/ur376715/Bgh/Bgh_seq/Genome/bgh_dh14_v4-2.bed

Genome_name=bgh_dh14_v4.fa


### promote conda environment 

source activate EDTA

cd $work_dir

cp $Genome_dir $work_dir

####################################

### run EDTA 
### --genome	[File]	The genome FASTA
### --cds	[File]	Provide a FASTA file containing the coding sequence (no introns, UTRs, nor TEs) of this genome or its close relative.
### --exclude	[File]	Exclude bed format regions from TE annotation. Default: undef. (--anno 1 required).
### --sensitive	[0|1]	Use RepeatModeler to identify remaining TEs (1) or not (0, default).
### 			This step is very slow and MAY help to recover some TEs.
### --anno	[0|1]	Perform (1) or not perform (0, default) whole-genome TE annotation after TE library construction.
### --evaluate	[0|1]	Evaluate (1) classification consistency of the TE annotation. (--anno 1 required). Default: 0.
###			This step is slow and does not affect the annotation result.
### --threads|-t	[int]	Number of theads to run this script (default: 4)
### --overwrite	[0|1]	If previous results are found, decide to overwrite (1, rerun) or not (0, default).

### the path of EDTA.pl has to be indicated. If the EDTA package were installed via conda, EDTA will locate at */conda_evn_name/share/EDTA/EDTA.pl
### this step is slow and will take around 24 hours

perl /home/ur376715/miniconda3/envs/EDTA/share/EDTA/EDTA.pl --genome $Genome_dir --cds $Genome_CDS_dir --exclude $Genome_BED_dir --overwrite 0 --sensitive 1 --anno 1 --evaluate 1 --threads 40 


####################################

### integrate repbase into repeatmasker

Repbase_dir=/hpcwork/ur376715/Data_base/Repbase_2018_Libraries/RMRBSeqs.embl

cp $Repbase_dir $RM_dir

### configure repeatmasker. This step will integrate repbase into repeatmasker database (combined with Dfam database)

cd $RM_dir

cd ..

### configure repeatmasker library

./configure

####################################

### assign names to Unknow consensus sensequence from EDTA library

### EDTA will generate a consensus library, which was used for annotation. The problem is that many consensus sequence don't have name

cd $work_dir
 
### -s  Slow search; 0-5% more sensitive, 2-3 times slower than default
### -pa(rallel) [number]
###        The number of sequence batch jobs [50kb minimum] to run in parallel.
###         RepeatMasker will fork off this number of parallel jobs, each
###         running the search engine specified. For each search engine
###        invocation ( where applicable ) a fixed the number of cores/threads
###         is used:

###           RMBlast     4 cores
###           ABBlast     4 cores
###           nhmmer      2 cores
###           crossmatch  1 core

###         To estimate the number of cores a RepeatMasker run will use simply
###         multiply the -pa value by the number of cores the particular search
###         engine will use.

RepeatMasker ${Genome_name}.mod.EDTA.TElib.fa -species fungi -s -pa 8 

### extract Unknow or unknow consensus ID

cat ${Genome_name}.mod.EDTA.TElib.fa.out | awk '/Unknown|unknown/' > ${Genome_name}.mod.EDTA.TElib.fa.out_Unknown_unknown.out

cat ${Genome_name}.mod.EDTA.TElib.fa.out_Unknown_unknown.out | awk 'BEGIN{FS=FS; OFS="\t"} {print $5,$11}' | awk 'BEGIN{FS="[# \t]"; OFS="\t"}{print $1"#"$2, $1"#"$3}'| uniq > ${Genome_name}.fa.mod.EDTA.TElib.fa.out_Unknown_unknown_ID_uniq

### uniq function can't totally eliminate duplicate ID. Then use uniq -d option to check if the first column is identical. If not, according to uniq file to delete duplication manually

cat *uniq | awk '{print $1}' | uniq -d > ${Genome_name}.fa.mod.EDTA.TElib.fa.out_Unknown_unknown_ID_uniq_checking_list

### according to the checking_list file to modify ${Genome_name}.fa.mod.EDTA.TElib.fa.out_Unknown_unknown_ID_uniq again 

### assign name to Unknon or unkown EDTA consensus library

awk '                          ##Starting awk program here.
FNR==NR{                       ##FNR==NR is condition which will be TRUE when 1st Input_file named ref.txt will be read.
  a[$1]=$2                     ##Creating an array named a whose index is $1 and value is $2 of current line.
  next                         ##next will skip all further statements from here.
}                              ##Closing BLOCK for FNR==NR condition here.
($2 in a) && /^>/{             ##Checking condition if $2 of current line is present in array a and starts with > then do following.
  print ">"a[$2]               ##Printing > and value of array a whose index is $2.
  next                         ##next will skip all further statements from here.
}
1                              ##Mentioning 1 will print the lines(those which are NOT starting with > in Input_file seq.fa)
' ${Genome_name}.fa.mod.EDTA.TElib.fa.out_Unknown_unknown_ID_uniq FS="[> ]" ${Genome_name}.mod.EDTA.TElib.fa > ${Genome_name}.mod.EDTA.TElib_after_repeatmasker_classification.fa


##########################################

### combine consensus sequence with fungi repbase sequence together

cat ${Genome_name}.mod.EDTA.TElib_after_repeatmasker_classification.fa Repbase_fungi.fa > ${Genome_name}.mod.EDTA.TElib_after_repeatmasker_classification_combined_with_fungi_repbase_consensus.fa

### run CD_HIT_EST to remove redundancies >80% similarity
### cd-hit can be installed by conda "conda install -c bioconda cd-hit"
### -c	sequence identity threshold, default 0.9 this is the default cd-hit's "global sequence identity" calculated as:
### 	number of identical amino acids or bases in alignment divided by the full length of the shorter sequence
### -M	memory limit (in MB) for the program, default 800; 0 for unlimitted;
### -T	number of threads, default 1; with 0, all CPUs will be used


cd-hit-est -i ${Genome_name}.mod.EDTA.TElib_after_repeatmasker_classification_combined_with_fungi_repbase_consensus.fa -o ${Genome_name}.mod.EDTA.TElib_after_repeatmasker_classification_combined_with_fungi_repbase_consensus_cd_hit_est_80.fa -c 0.8 -n 5 -M 10000 -T 0

### run repeatmasker based on the combined consensus library 

### -html Creates an additional output file in xhtml format.

### -gff  Creates an additional Gene Feature Finding format output

### -xm Creates an additional output file in cross_match format (for parsing)

### -poly  Reports simple repeats that may be polymorphic (in file.poly)


RepeatMasker $Genome_dir -lib ${Genome_name}.mod.EDTA.TElib_after_repeatmasker_classification_combined_with_fungi_repbase_consensus_cd_hit_est_80.fa -pa 30 -s -a -gff -html -poly -xm

### parseRM to check the result
### -v,--v (BOOL)
###         Verbose mode, make the script talks to you
###         Print the version if only option

### -l,--land (STRING)
###         To generate additional outputs that can be used to make 
###         landscape graphs by repeat, by family and by class (3 files).
###         Two values should be given, separated by a comma: <max>,<bin>
###            <max> is the value of the last bin, where everything 
###                  higher than this %div or My will be placed
###            <bin> is the size of the bins
###         If -m is set, then numbers here HAVE TO BE in My instead of %divergence
###         Typically: -l 50,1 (if %div are corrected for CpGs then values tend to be higher)
###         See the end of this help for more examples.

### -i,--in (STRING) 
###         RepeatMasker output .out or RepeatMasker output .align 
###         (use -d and provide a directory path here to parse several files)
###         Use of .align is best for -l, because the %div corrected 
###         for higher mutation rate at CpG sites and the Kimura 2-Parameter 
###         divergence metric is in the .align files (thus, the graphs are smoother)
###         However, the script will treat separately repeats that overlap 
###         (merged in the .out, with one name kept): this means that for -p, 
###         the amount of DNA masked by each of the repeats will be much higher 
###         than if a .out file is parsed.



perl /hpcwork/ur376715/Softwares/Parsing-RepeatMasker-Outputs-master/parseRM.pl -i ${Genome_name}.align -l 50,1 -v

################################## note
 
### the tools used for extracting species specific sequences "queryRepeatDatabase.pl" from Repbase was not available for repeatmasker whose version is greater than 4.1.1-0 
### for this reason, the lower version of repeatmasker was downlaoded and used for repbase filter. 

conda install -c bioconda repeatmasker=4.1.0-0

### copy repbase into this repeatmasker

### the species specific sequences were extracted 

ur376715@login18-1:~/miniconda3/envs/TEtranscript/share/RepeatMasker/util[685]$ perl queryRepeatDatabase.pl -species fungi > /work/ur376715/TE_analysis/Repbase_fungi.fa

ur376715@login18-1:~/miniconda3/envs/TEtranscript/share/RepeatMasker/util[686]$ perl queryRepeatDatabase.pl -species ascomycota > /work/ur376715/TE_analysis/Repbase_ascomycota.fa


##################################














