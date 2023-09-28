# NOTES ON HYBRIDIZATION ANALYSES
## Bayescenv Analysis

1.	Convert the structure file populations_n5nb_str3.txt (Wren_1/Clines) to Bayescenv format using PGDSpider. This was the file for a second run of data in Structure using sampling locations as priors. Code for getting the sampling location codes is in the Clines Rscript. I replaced in Excel zeros with -99 to help PDGSpider to identify missing data.

2.	Get the environmental variable from Chelsa in Clines Rscript, and write the scale vector as a text file.

3.	Once I downloaded the compressed files from https://github.com/devillemereuil/bayescenv, I copied these files into my Linux folder. I changed the permission rights before decompress (see command lines for that). I decompressed the folder (see command lines for that). I copied and pasted the data files into the decompressed folder ./linux/bayescenv-1.1/bin/linux64/. Finally, I change the permission rights of the data files.

sudo chmod -R 777 ./AMPp.txt

4.	I set the priors for pi (probability of departure from the neutral model) at 0.1, and the prior for p(probability that the selection is driven by the given environmental variable) at 0.5. According to the author these are conservative settings. I used the default priors for g and alpha.

5.	Run the software in Ubuntu with the following command line. I sued the default MCMC setting of the software.
./bayescenv bayescenv.data -snp -env AMPp.txt -threads 10 -o run1 -pr_jump 0.1 -pr_pref 0.5
6.	I had some issues running the environment file in Ubuntu. It worked after I modified slightly the file in Ubuntu using the nano command and saving it. I guess there is an issue with Windows/Linux compatibility. Check the EOL Conversion options in Notepad++ (check the instructions for this).
