# ICEBERG
Scripts to execute the ICEBERG pipeline as described in https://doi.org/10.1093/nar/gkae180

This is a fully executable version of ICEBERG, you just put it in a folder with your samples and run it. It is set up as an executable bash script, that should run you through from the raw files to the final ICEBERG peaks and files for visualization, including basic heatmaps, bigwigs, etc. 

Your samples should all be in the same folder with the script, preferably as raw fastq.gz files, named as Target_1, Target_2, etc. and you should have an already processed bam file for your IgG/negative control for peak calling. 

You should have one environment with all of the programs needed to run the script, except for deeptools. I have provided a file that can be used to create this environment with the proper programs (package-list-cr.txt). Create a second environment for deeptools (version >= 3.5.4), this is due to incompatabilities with other programs in the main environment. You must provide the names of these two environments in the beginning of the script. 

You will need to go in and edit the variables at the beginning to reflect your base naming structure, number of parallel threads to use, number of replicates, and path to the files and genomes, and there is also an option to change the MACS2 q-value threshold if needed. You will need an appropriate blacklist/suspect list for artifact region filtering, and a header sam file (header.sam) - these are included in this repository for hg38 and mm10 genomes.  

This script uses our usual Cantù lab in-house settings and programs from trimming, mapping, etc. I set it up to use MACS2 like we do in the paper, but that can also be changed to whatever peak caller you like. The script will automatically calculate the replicate with the lowest mapped reads, and then the down sampling will happen according to that number. It will create 3 aggregates, and the ICEBERG peaks will be defined as the union (2 of 3) between them. It will also call peaks after adding each subsequent replicate (and add replicates in a random order each time), and log this in a text file plus make a simple line graph to mimic the build curves in the paper. It makes a .bw file for visualization that is the average of the 3 aggregates, and also plots a basic heatmap for the Iceberg peaks in all the replicates and in the final .bw. 
