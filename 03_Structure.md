# NOTES ON HYBRIDIZATION ANALYSES
## Population Structure with STRUCTURE

Besides ADMIXTURE I will also work with STRUCTURE, which uses a Bayesian framework for estimating the probability of a sample to belong to a group. The first problem I found working with STRUCTURE is that I must manually made some minor editing in a text editor to the data file. STACK can export the data as STRUCTURE format, but it includes some information of data and time in the first line. STRUCTURE cannot process this, so we must delete the first line manually in a text editor. However, when I save it in windows as a text file, it does it in a format incompatible with Linux and not possible to run in the cluster. So, the alternatives are run the data in the local computer which might take a long time or even be unfeasible. Alternatively, it is possible to save as UNIX OS in Notepad ++ going to Edit > EOL Conversion > UNIX/OSX Format. The other solution might be converting the format to a Linux text file in the terminal with the following lines.
module load perl
perl -p -e s/\r$// < populations_n5nb.txt > populations_n5nbx.txt
We can also convert the .vcf or the .ped file from STACK to a STRUCTURE format using PGDSpider. http://www.cmpg.unibe.ch/software/PGDSpider/#Download_and_Installation_Instructions

STRUCTURE requires three files. The input file with the genomic data, mainparams file with the number of samples, loci, burn-in, replicates and other file-format parameters, and the extraparams file with info for the priors.
Input File. It is always important to check the format of the input file we get. The first row contains the marker names. Each individual is represented by two rows; if the individual is homozygous at a marker then the two rows will be identical, if the individual is heterozygous at a marker then the two rows will differ for that marker, if the individual is haploid (male) then the second row of data for that individual will consist of missing data at each marker 
1 represents A 
2 represents T 
3 represents C 
4 represents G 
-9 (zero or any number not present in the genomic data) represents missing data (due to an individual being haploid, or due to no genotype recorded for that individual at that marker) 
Main and Extra Parameter Files. Besides the main data file, we will need a mainparams and extraparams file which are files without extensions that contain the parameter for the model. Info for mainparams and extraparams
http://cbsuapps.tc.cornell.edu/extraparams.txt
http://cbsuapps.tc.cornell.edu/mainparams.txt
http://cbsuapps.tc.cornell.edu/blast_s_helpwin.aspx?hid=mainparams&app=structure
http://cbsuapps.tc.cornell.edu/blast_s_helpwin.aspx?hid=extraparams&app=structure
https://github.com/sanger-pathogens/structure/blob/master/mainparams
1.	Once we get the data input and the parameters files ready, we move them into the cluster and execute the following bash script to run STRUCTURE from K = 1 to K = 7 and 30 replicates. We can also just modify and use the script Wren_14_Struc_n5nb.sh. I choose to set the burn-in at 100000 and the interactions at 1000000 following the recommendation of Gilbert et al. (2012).

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3500mb              
#SBATCH --time=4-00:00:00   
#SBATCH --job-name=Wren_n5nb 
#SBATCH --qos=robinson-b   #Use low performance no used resources
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=Wren_n5nb_%A_%a.out      
#SBATCH --array=9-30
pwd; hostname; date


module load structure/2.3.3

### Run the loop of runs for this task.
for ((k=1; k<=7; k++)); do
 
structure -K $k \
-i /blue/robinson/ldmontalvo/wrens_1/Struc/populations_n5nbx.txt \
-o /blue/robinson/ldmontalvo/wrens_1/Struc/n5nb_2/wren_struc_n5nb_k${k}_r${SLURM_ARRAY_TASK_ID} \
-m /blue/robinson/ldmontalvo/wrens_1/Struc/mainparams \
-e /blue/robinson/ldmontalvo/wrens_1/Struc/extraparams
done

date
Structure documentation: 
https://www.ccg.unam.mx/~vinuesa/tlem09/docs/structure_doc.pdf
https://web.stanford.edu/group/pritchardlab/software/readme/readme.html

2.	 An excellent GitHub site that explain the process (including parameters) to use STRUCTURE  in a cluster is at https://github.com/mossmatters/Structure-Pipeline.
3.	In the case of my data, I also run a model with LOCPRIOR =1 and LOCISPOP=1 option. This is advised for data that dont have strong structure. This option uses the preassigned geographical population as prior for the posterior assignment. We also keep the NOADMIX=0 to allow for admixture in the assignment. We use the correlated frequencies model since we assumed common ancestries for most the individuals (FREQSCORR=1) (Porras-Hurtado et al., 2013; Pritchard et al., 2003).

## Finding K with Harvester.
4.	Once we get the results from STRUCTURE we can download the python scripts for Harvester https://github.com/dentearl/structureHarvester and run it. Alternatively, we have to zip-compress the files form structure and upload them into the website.
5.	The Delta k method detects the uppermost level of population structure when several hierarchical levels exist. You can try to perform a hierarchical analysis of population genetic structure as described in Coulon et al 2008 (Molecular Ecology). The idea is to perform supplementary STRUCTURE analyses for each of your Ks independently, and to repeat the process until you find values of k=1. This may allow you to characterize substructure in your data. We must consider if we have enough sample size to partition the data into smaller hierarchical subset of the data.

python structureHarvester.py --dir=/blue/robinson/ldmontalvo/wrens_1/Struc/test --out=/blue/robinson/ldmontalvo/wrens_1/Struc/test --evanno --clumpp

6.	Running Harvester using Windows cmd. To run the harvester in Windows, first download the structureHarvester.py and putting in the directory where you want to work. Put the folder with all the outcome files from Structure (ex. In my case n5nb_2) in the same directory where the structureHarvester.py file is. Create a folder for the outcomes (ex. Outcome1).

cd dir # Directory where the structureHarvester.py, the folder with Structure outcome and empty folder for results are.

structureHarvester.py --dir=.\n5nb_2 --out=./Outcome.k1to10_r30 --evanno --clumpp

Harvester creates files called K*.indfile consisting in the proportions assigned by structure for Kn for all the replicates or runs in a single file that we can use in CLUMPP.

7.	It turns out that K=7 had the highest delta K. this means the upper level of clustering in the data is seven groups. However, the barplot for k=2 for one of the replicates suggested there is Isolation by Distance based on the gradual change of genetic composition (Fig. 7). Under this scenario, Structure do not work well detecting subgroups or groups since it is hard to establish limits of the groups. We can try one of the spatial explicit software that use geographic coordinates such as:



![Presentation1](https://github.com/ldmontalvo/Landscape-Genomic-Wrens/assets/67880755/15ef10d1-21b3-4c2e-8acf-8fb37c31ee36)











8.	Some options for incorporating spatial info:
TESS: http://membres-timc.imag.fr/Olivier.Francois/
BAPS: http://www.helsinki.fi/bsg/software/BAPS/
Geneland: https://i-pri.org/special/Biostatistics/Software/Geneland/

## Summarizing Q matrices with CLUMPP.
9.	Harvester gives us an indfile which put together all the runs or replicates for each values of K. We can use this indfile from Harvester as input for CLUMPP.

CLUMPP K1.paramfile -k 1
