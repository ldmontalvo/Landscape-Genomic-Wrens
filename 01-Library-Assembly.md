# NOTES ON HYBRIDIZATION ANALYSES
## Library construction and De Novo Assembly

## Trim and Filtration of Low Quality Sequences Using fastx_toolkit module
1.	Use the bash file Wren_3_trim (or paste the script below) to cut the last bases in the sequences that had more error, base on the Phred scores (first graphic in the fastqc results) you can see how much error did you get at the ends of the sequences. The code here is set to cut at 147 bases. This means we let out the last three bases of the sequences (150 bp long). We got very good reading at the end.

#!/bin/bash
#SBATCH --ntasks=1                        
#SBATCH --mem=1gb                          
#SBATCH --time=07:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load fastx_toolkit

echo "Trim 3' end of all reads to a length of 147 bp (FASTX Trimmer)"

fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112439_HW3FFAFXY_DM1_ATCACG_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM1_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112440_HW3FFAFXY_DM2_CGATGT_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM2_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112441_HW3FFAFXY_DM3_TTAGGC_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM3_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112442_HW3FFAFXY_DM4_TGACCA_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM4_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112443_HW3FFAFXY_DM5_ACAGTG_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM5_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112444_HW3FFAFXY_DM6_GCCAAT_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM6_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112445_HW3FFAFXY_DM7_CAGATC_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM7_t.fastq &
fastx_trimmer -f 1 -l 147 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/11423_3270_112446_HW3FFAFXY_DM8_ACTTGA_R1.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM8_t.fastq

nohup bash Wren_3_trim.sh &

date

2.	Use the bash file Wren_4_filter (or paste the script below) to eliminate some sequences that had low quality scores. We did these twice. One cut any sequences below 10 and then, using these cut files, we eliminate again the ones below 20. For this task I added cpus-per-task=8 and --mem-per-cpu=1gb to speed up the time needed. These parameters mean that I will use eight cpus for these tasks and each one will use 1 gb of ram to process, being 8 times faster than what was specified in the previous script. By using the parameters of the previous script (in Wren_3_trim), this process takes up to seven hours. With the new parameters, the process takes 28 mins. For an explanation of Phred quality score see: https://en.wikipedia.org/wiki/Phred_quality_score

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=02:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load fastx_toolkit

echo "Filter sequences with less than 10 and 20 Phred scores"


##Eliminate sequences where there is a sinlge Phred score below 10 and then sequences where 5% of reads have a with Phred quality scores below 20

label of "tff" indicates that the second round of quality filtering is completed

#100% of the bases in a sequence must have a score of higher than 10 for the sequence to be kept
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM1_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM1_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM2_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM2_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM3_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM3_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM4_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM4_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM5_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM5_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM6_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM6_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM7_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM7_tf.fastq &
fastq_quality_filter -q 10 -p 100 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM8_t.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM8_tf.fastq

#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM1_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM1_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM2_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM2_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM3_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM3_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM4_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM4_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM5_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM5_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM6_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM6_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM7_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM7_tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/AS_LDM8_tf.fastq -o /blue/robinson/ldmontalvo/wrens_1/AS_LDM8_tff.fastq


nohup bash Wren_4_filter.sh &

date

## De Novo Assembly with STACK
3.	Next step is the demultiplexing. Use the script Wren_5_demultiplex.sh or copy and paste the code below. The indexes used can be seen below. They are shown at the end of the Illumina file names. With the parameters cpus-per-task=8 and --mem-per-cpu=1gb, the process took 11 min. It is important to keep a consistent format of codes or names for the samples, be aware of spaces that might not show up and the end of a line. Try to keep the names as consistent as possible. You will change the part highlighted in the DNA string in the code with the correspond index illumina code. Data from illumina will bring the code (the six bases code per index) in the name of the folder (e.g 11423_3270_112439_HW3FFAFXY_DM1_ATCACG_R1.fastq.gz). This index is standardized and are frequently the same in the same order. The R1 in the name of the file means the data is just one direction of read. When the sequence is done in both directions, there will be two compressed folders, instead of one, one folder ending in R1 and the other in R2.

Indexes used (illumina)
ATCACG
CGATGT
TTAGGC
TGACCA
ACAGTG
GCCAAT
CAGATC
ACTTGA

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=10:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date


mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM2raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM3raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM4raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM5raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM6raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM7raw
mkdir /blue/robinson/ldmontalvo/wrens_1/AS_LDM8raw


mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM1_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM1_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM2_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM2_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM3_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM3_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM4_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM4_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM5_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM5_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM6_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM6_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM7_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM7_tff.fastq
mv /blue/robinson/ldmontalvo/wrens_1/AS_LDM8_tff.fastq /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw/AS_LDM8_tff.fastq

mkdir /blue/robinson/ldmontalvo/wrens_1/demultfilter

module load stacks


process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM1raw -b /blue/robinson/ldmontalvo/wrens_1/index1.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM2raw -b /blue/robinson/ldmontalvo/wrens_1/index2.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM3raw -b /blue/robinson/ldmontalvo/wrens_1/index3.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM4raw -b /blue/robinson/ldmontalvo/wrens_1/index4.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM5raw -b /blue/robinson/ldmontalvo/wrens_1/index5.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM6raw -b /blue/robinson/ldmontalvo/wrens_1/index6.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM7raw -b /blue/robinson/ldmontalvo/wrens_1/index7.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
process_radtags -p /blue/robinson/ldmontalvo/wrens_1/AS_LDM8raw -b /blue/robinson/ldmontalvo/wrens_1/index8.txt -o /blue/robinson/ldmontalvo/wrens_1/demultfilter -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina


#move the demutliplex_AS_LDM.sh and index files to the working directory using fetch.

nohup Wren_5_demultiplex.sh &

date

4.	We will select the samples with enough number of read (>100000). We can see the number of reads in the process_radtag.ASLDM1raw.log file. There will be a file for every index. You can also run the following command in the terminal with the appropriated directory: wc /blue/robinson/ldmontalvo/wrens_1/demultfilter/*.fq -l (Fig 3). Copy and paste this information in Excel and divide the column with the number of read by 4 (4 reads). This will be the number of reads we will consider to select the samples (Fig. 4)

5.	Select around 12 of the samples with the highest number of reads and make a text file tab delimited with the names of the samples. We will run de novo assembly several times to optimize the parameters M and n. -M is basically the number of mismatches allowed between stacks within individuals and -n is the number of mismatches allowed between stacks between individuals. We try different parameters for -M and -n (from 1 to 9). Change the corresponding value for n so that we follow the n=M rule. Fix M = n, from 1 to 9 (the number of mismatches between 2 alleles in either a heterozygote M, or in population n) and keep m = 3 (stack depth/number of identical reads required to initiate a new allele). We are testing all these parameters with a minimum percentage of individuals in a population required to process a locus for that population (-r 0.8). It’s been showed that at least 80% of the population should present a specific locus to be included, known as the 80% rule or r80 loci (Paris et al., 2017). We can use the script Wren_6_denovo_test.sh or copy and paste the code below. For a more detailed information about de novo assembly with stacks: http://catchenlab.life.illinois.edu/stacks/manual/#pfiles
6.	I started the analysis with 158 samples, 20 of them were repeated. These were group 5 and 8 that were prepared twice because of the low concentration of DNA. I really had 138 samples. From the 158 samples, I took off 16 with low readings. These included 2 samples already in the repeated group. So, I took off 14 samples. I additionally took off 3 samples of C bruneicapillus that I suspected had some cross-contamination. These summed up 17 samples. So, I ended up with 121 samples for the Denovo assembly. 
Removed samples
Samples of C. brunneicapillus with low reads B008, B009
Other samples with low reads: 74, 33, 78, 84, 87B, 87, 91, 81, F018, 90, 14, 89, F027, 70
Samples of C. brunneicapillus with cross contamination or mislabeling: B004, B010, B011

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=10:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load stacks

##Run 1

mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run1
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 1 -n 1 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run1
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run1 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8

#note -T or -t = # of threads/CPUs to use
#for populations: -P = path to directory containing the Stacks files (output from the denovo run)
#-M = path to the population map (in thsi case same as used in the denovo trial)
#-r = minimum percentage of individuals in a population req to have a locus to process the locus

##Run 2
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run2
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 2 -n 2 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run2
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run2 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8


##Run 3
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run3
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 3 -n 3 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run3
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run3 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8


##Run 4
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run4
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 4 -n 4 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run4
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run4 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8


##Run 5
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run5
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 5 -n 5 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run5
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run5 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8


##Run 6
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run6
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 6 -n 6 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run6
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run6 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8


##Run 7
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run7
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 7 -n 7 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run7
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run7 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8



##Run 8
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run8
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 8 -n 8 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run8
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run8 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8


##Run 9
mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_run9
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -m 3 -M 9 -n 9 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_run9
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_run9 -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_test.txt -t 15 -r 0.8



nohup Wren_6_denovo_test.sh &


date

7.	We will get a graph of number of loci per each parameter. We will use the script Wren_7_getp_test.sh or the code below. This code will generate nine different files (one per parameter) in your directory with the number of SNPs and loci per SNP.

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=00:20:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load stacks

cd /blue/robinson/ldmontalvo/wrens_1/

for i in 1 2 3 4 5 6 7 8 9
do
  stacks-dist-extract denovo_run${i}/populations.log.distribs snps_per_loc_postfilters > ${i}_snp_distribution.tsv
done

nohup Wren_7_getp_test.sh &

date

8.	We will switch for a while to R and use the script Par_Selec or the code below to make a graph of number of SNPs per parameter (source: https://www.r-bloggers.com/parameter-testing-in-stacks-snp-extraction-and-visualization-in-r/).

count <- 1
files <- list.files(path = "C:/Users/Owner/Desktop/Dropbox/Thesis/Molecular_Wrens/Wren_1/FastQC_Stacks", pattern="*_snp_distribution.tsv", full.names = T)
for (i in files[1:9]){
  table <- read.delim(i, skip=1, header=T)
  table$n_loci_percent<- table$n_loci/sum(table$n_loci)
  table$m<- count
  write.table(table, "distributions.tsv", append=T, row.names=F, col.names = F)
  snp_count <- data.frame("m"= count, "n_snps"=sum(table$n_loci))
  write.table(snp_count, "total_count.tsv", append=T, row.names=F, col.names = F)
  count <- count + 1
}

#total_count.tsv is used to display total SNP by parameter.

#Figure 5
library(ggplot2)
snp_count<-read.delim("total_count.tsv", sep=" ", header=F)	
names(snp_count)<-c("m", "n_snps")
snp_count$m<-as.factor(snp_count$m)
p <- ggplot(data=snp_count, aes(x=m, y=n_snps))
p + geom_point(size=4)

#Figure 6
snp_table<-read.delim("distributions.tsv", sep=" ", header=F)
names(snp_table)<- c("n_snps","n_loci", "n_loci_percent", "m") 
snp_table$n_loci_percent<-snp_table$n_loci_percent*100
snp_table$n_snps<-ifelse(snp_table$n_snps < 9, snp_table$n_snps, "9 +")
snp_table$n_snps<-as.factor(snp_table$n_snps)
snp_table$m<-as.factor(snp_table$m)

q<-ggplot(data = snp_table) + 
  geom_col(aes(x=n_snps, y=n_loci_percent, fill=m), position="dodge") + theme_classic()
q

9.	We can see after m=4 or m=5 we really dont get much more SNPs. So, we can use 5 as M and n parameters for the whole set of samples.





![Slide1](https://github.com/ldmontalvo/Landscape-Genomic-Wrens/assets/67880755/16b659cd-fb98-4f87-8720-642fbeff06cb)
![Slide2](https://github.com/ldmontalvo/Landscape-Genomic-Wrens/assets/67880755/04a7e7eb-a6fd-42b2-b2aa-a2df52deb35e)









10.	Now we have our M and n number and know that over M=n=5 we dont get much more SNPs. In other words, if we increase the mismatches or error to accept a SNP over five, we dont get necessary much more SNPs. We also need a good coverage of the loci. By convention, we need 20x more genomic material, so to get >95% probability of correctly calling a SNP, one needed at least about 20x genomic coverage. Some studies showed that to get an approximate of 20X of coverage we need at least 100,000 reads (Hoffberg et al., 2016). We will use samples that have at least this number of reads and make a text file as in step 17 with these samples. We will repeat the step 21, but with the parameter M=n=5 for all the samples with more than 100,000 reads. In this step, we can filter out any sample we do not want by not including them in the text file. For example, we can leave out one species or population. We run script Wren_8_denovo_final.sh or the code below. With the parameters in the code, it took. We change some output parameters in the populations command to get vcf file for further analyses in R. For more information of the parameters for SNPs calling filtering see: https://catchenlab.life.illinois.edu/stacks/comp/populations.php.

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=10:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load stacks

mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_final

denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_final.txt -m 3 -M 5 -n 5 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_final
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_final -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_final.txt -r 0.8 -p 1 -t 15 --min_maf 0.05 --write_single_snp --structure --vcf --

nohup Wren_8_denovo_final.sh &

date
11.	We must find the best combination of parameters to filter bad quality of sequences. We will run two assembles, one dividing the samples in two populations (changing the popmap file) and another excluding the outgroup C. brunneicapillus for the IBE analyses. Since we are using two populations, we use other parameters here. The setting populations provides filtering options to only include loci that occur at certain frequencies in each population. If you set the SNP filtering to two or more populations, we can use the parameters -R (minimum percentage of individuals across populations required to process a locus) and --max-obs-het (specify a maximum observed heterozygosity required to process a nucleotide site at a locus). If we have only samples of haploid cells, we set up the maximum accepted heterozygosity at 0. We set -R 0.5 if we considered all individual a single population and eliminated loci that were missing from more than 50% of individuals. We can use the file Wren_10_denovo_n5_2pop.sh or use the code below. Notice that in the command population we change the parameter -p 2 and add and add the parameters -R 0.5 and --max-obs-het 0.7.
According to Jim Austin, having r = 0.8 might be a little too conservative. This is 80% of your samples within a population (however you’re defining one) have to have the locus sequenced. The parameter p=1 might be not useful. It is essentially allowing loci in that are only found in one population. So, they cant be compared to anything. Alternatives are to adjust the popmap file to be either all birds in one pop then run through populations step, or alternatively 2 populations one in-group species and one for brunneicapillus. This way, it will be easier to separate out by species or pop/location downstream but it will reduce the chance of introducing bias in the STACKS step. (for example the clear outlier status of C. brunn in the PCA - If all loci in only 1 population (p=1) then there might be a bunch of loci only in that species, or missing from that species, then the PCA could be driven by that.

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=3gb              
#SBATCH --time=05:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load stacks

mkdir /blue/robinson/ldmontalvo/wrens_1/denovo_n5_2pop
denovo_map.pl --samples /blue/robinson/ldmontalvo/wrens_1/demultfilter --popmap /blue/robinson/ldmontalvo/wrens_1/popmap_wren_n5_2pop.txt -m 3 -M 5 -n 5 -T 15 -o /blue/robinson/ldmontalvo/wrens_1/denovo_n5_2pop
populations -P /blue/robinson/ldmontalvo/wrens_1/denovo_n5_2pop -M /blue/robinson/ldmontalvo/wrens_1/popmap_wren_n5_2pop.txt -r 0.8 -R 0.5 -p 2 -t 15 --min_maf 0.05 --max-obs-het 0.7 --write_single_snp --structure --vcf --genepop --hzar --fasta-loci --fasta-samples --plink

nohup Wren_10_denovo_n5_2pop.sh &

date
12.	The link below explained the difference between VCF file (variant sites only) and .fa files (invariant and variant sites):
https://groups.google.com/g/stacks-users/c/BJxvnQ79OG0

Hardy-Weinberg.
13.	 Vcftools is probably one of the commonly used tools for pruning genotypes data for linkage disequilibrium and Hardy-Weinberg disequilibrium. First, Install vcftools in Ubuntu:
sudo apt-get install -y vcftools
14.	We use the vcf file we got from Stack and the procedure and script from http://www.ddocent.com/filtering/. First we install vcflib with:
conda install -c bioconda vcflib
15.	We convert to SNPs:
vcfallelicprimitives populationsn5snps.vcf --keep-info --keep-geno > populationsn5recode.vcf
vcftools --vcf populationsn5recode.vcf --remove-indels --recode --recode-INFO-all --out forhwe
16.	We use the perl script by Chris Hollenbeck at http://www.ddocent.com/filtering/. We used a threshold of p-value of 0.05 and filtered out SNPs that are below that p-value threshold in 0.5 proportion of all populations. In other words, we excluded any locus that was out of HWE (a = 0.05) at >50% of these sites.
./filter_hwe_by_pop.pl -v forhwe.recode.vcf -p poplab.txt -o vcftools.hwefilter -h 0.001 -c 0.5
17.	We can use this file for later analysis.
Linkage Disequilibrium with vcftools
These steps were run in the cluster directly in the terminal. This can be also run in my local computer, but it will need the installations of the applications. We will try two ways to prune the vcf file from STACK to get rid of the loci or SNPs that are under linked.
Linkage Disequilibrium Pruning with vcftools
18.	Run the following lines in the terminal or directly in the local computer if the applications are installed.
module load vcftools
vcftools --vcf pop.snps.vcf --hap-r2 --min-r2 .7 --ld-window-bp 50000 --out stats_ld
19.	I used vcftools to apply other filters. I could use it to filter minor allele frequency (maf) which for my data I already did in Stack in the populations command.
vcftools --vcf pop.snps.vcf --maf .3 --recode --out stats_maf
20.	It is possible also to prune the data for alleles that are not in HWE with the following line.
vcftools --vcf pop.snps.vcf --hwe 0.01 --recode --out n5nb_filter
This line takes out any allele that diverge from HWE or is significantly different from HWE expectation with a p-value of 0.01. With this criterion I got 2879 SNPs.
Other options for filtering:
--min-r2 .7 is the threshold for the minimum correlation coefficient to discard SNPs.
--ld-window-bp <integer> This optional parameter defines the maximum number of physical bases between the SNPs being tested for LD in the "--hap-r2", "--geno-r2", and "--geno-chisq" functions.
In this data set there wasn’t any correlation coefficient over 0.7. 
--chr <chromosome>
--not-chr <chromosome>
Includes or excludes sites with identifiers matching <chromosome>. These options may be used multiple times to include or exclude more than one chromosome.
--geno-chisq Outputs a file reporting the r2. This options is used with multiallelic data (more than to alleles).
--min-r2 <float> This optional parameter sets a minimum value for r2, below which the LD statistic is not reported by the "--hap-r2", "--geno-r2", and "--geno-chisq" functions.
More information of VCFTools at http://vcftools.sourceforge.net/man_latest.html
Descriptions of the options https://vcftools.github.io/man_0112b.html
A good guide for VCFtools is http://www.ddocent.com/filtering/
21.	Additionally, we can use vcftools to subset the vcf file. We need a txt file with one columns no headers with ID of samples. One ID per row.
vcftools --vcf n5nb_filter.vcf --remove CalcetaSamples.txt --recode --recode-INFO-all --out n5nb_ncal_fil
Linkage Disequilibrium with plink
22.	Filter Out SNPs to Remove Linkage Disequilibrium (LD): SNPs in high LD with each other contain redundant information. More worrisome is the potential for some regions of the genome to have a disproportionate influence on the results and thus distort the representation of genome-wide structure (Liu et al., 2020). A standard approach to address this issue is to filter out SNPs based on pairwise LD to produce a reduced set of more independent markers. Here we use plinks commands to produce a new LD-pruned dataset (Liu et al., 2020). We use the plink.ped and plink.map files from the last step in Stack. The filtering is done in two steps. First, we set the parameters to prune the database and get the list of SNPs that will be pruned. Then, we export the pruned data. Here, I used the command –file for .ped and .map files other wise use –bfile for .bed files formats as inputs. Check the table with commands for plink for more options (http://zzz.bwh.harvard.edu/plink/reference.shtml ).
For more details about the command Plink see  http://zzz.bwh.harvard.edu/plink/summary.shtml#prune or  https://www.cog-genomics.org/plink/1.9/ld.

We can run the bash script below (the Wren_12_LD_n5.sh file). The analysis with the parameters for r2 = 0.2 did not find any SNP pair over that threshold. The analysis took ~2:10 min, so the parameters for hipergator memory were underused.  So, we keep the same original .ped and .map file for the next step in ADMIXTURE. We just have to change the directory or name of the .map file to run any other assembly. 

#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=3gb              
#SBATCH --time=05:00:00   
#SBATCH --job-name=check_test              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test_%j.log         
pwd; hostname; date

module load   plink/1.07

plink --file /blue/robinson/ldmontalvo/wrens_1/LD/populations --indep-pairwise 50 5 0.2
plink --file /blue/robinson/ldmontalvo/wrens_1/LD/populations --extract plink.prune.in --make-bed --out wren_LD_n5_2pop


date
Both vcftools and plink show no SNPs under linkage disequilibrium so I kept the same vcf file for the rest of analyses.
