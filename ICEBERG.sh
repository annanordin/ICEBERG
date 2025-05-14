#!/bin/bash
#
#SBATCH -J ICEBERG
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH --mail-user=
#SBATCH --mail-type=end

# General pipeline for CUT&RUN samples followed by ICEBERG pipeline. 
# Required programs: bbmap, bowtie2, samtools, bedtools, MACS2, deeptools. I have 2 different envs, separate for deeptools
# Required files: bowtie2 human genome, hg38/mm10 suspect list to remove artifacts, sam file containing only header information (header.sam), IgG.bam file for peak calling.

# Set number of threads
THREADS=

# Set number of replicates
REPLICATES=

# Set MACS2 q-value/p-value threshold (default is "-q 5e-2")
MACS2_THRESH="-q 5e-2"

# Set path to genome
GENOME=/path/to/genome

# Set path to suspect list
SUSPECTLIST=SuspectList_hg38.bed

# Set path to IgG bam file for peak calling
IgG=IgG.bam

# Set env variables
CUTNRUN=cutnrun
DEEPTOOLS=deeptools

set -e

# Check if required files exist
if [ ! -f "$SUSPECTLIST" ]; then
    echo "Error: Suspect list file $SUSPECTLIST not found!"
    exit 1
fi

if [ ! -f "$IgG" ]; then
    echo "Error: IgG file $IgG not found!"
    exit 1
fi

if [ ! -f "header.sam" ]; then
    echo "Error: Header file header.sam not found!"
    exit 1
fi

# Create directories for aggregates if they do not exist
mkdir -p Aggregate1 Aggregate2 Aggregate3

# Create a log file to record progress
LOGFILE="ICEBERG_log.txt"
echo "Starting pipeline execution: $(date)" > $LOGFILE

# Loop through samples
for sample in *_R1.fastq.gz; do

    source activate $CUTNRUN

    base_name=$(basename $sample _R1.fastq.gz)  # Extract base name of sample

    echo "Processing sample: $base_name" >> $LOGFILE

    # Trimming reads (bbduk)
        echo "Trimming for sample $base_name..." >> $LOGFILE
        bash bbduk.sh in=${base_name}_R1.fastq.gz in2=${base_name}_R2.fastq.gz \
            out=${base_name}_R1.fastq out2=${base_name}_R2.fastq \
            ref=adapters,artifacts literal=TATATATATATATATATATATATATATATATATATA,ATATATATATATATATATATATATATATATATATAT,GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC >> $LOGFILE 2>&1
        if [ $? -ne 0 ]; then
            echo "Error in bbduk for $base_name" >> $LOGFILE
            exit 1
        fi

    # Mapping reads to the genome (bowtie2)
    echo "Mapping reads to genome for sample $base_name..." >> $LOGFILE
    bowtie2 -p $THREADS --local --very-sensitive-local --no-unal --no-mixed --no-discordant --dovetail \
        -I 0 -X 500 -x $GENOME -1 ${base_name}_R1.fastq \
        -2 ${base_name}_R2.fastq \
        -S ${base_name}.bowtie2.sam >> $LOGFILE 2>&1
    if [ $? -ne 0 ]; then
        echo "Error in bowtie2 mapping for $base_name!" >> $LOGFILE
        exit 1
    fi

    # Converting SAM to BAM, sorting, and removing duplicates
    echo "Processing SAM to BAM for sample $base_name..." >> $LOGFILE
    samtools view -@ $THREADS -b ${base_name}.bowtie2.sam > ${base_name}.bowtie2.bam
    samtools fixmate -@ $THREADS -m ${base_name}.bowtie2.bam ${base_name}.bowtie2.fixmate.bam
    samtools sort -@ $THREADS ${base_name}.bowtie2.fixmate.bam -o ${base_name}.bowtie2.positionsort.bam 
    samtools markdup -r -@ $THREADS ${base_name}.bowtie2.positionsort.bam ${base_name}.bowtie2.markdup.bam 
    samtools sort -@ $THREADS ${base_name}.bowtie2.markdup.bam -o ${base_name}.bowtie2.final.bam 
    if [ $? -ne 0 ]; then
        echo "Error in SAM to BAM conversion for $base_name!" >> $LOGFILE
        exit 1
    fi

    # Removing suspect list overlaps
    bedtools intersect -abam ${base_name}.bowtie2.final.bam -b $SUSPECTLIST -v | samtools sort -o ${base_name}.bam 
    if [ $? -ne 0 ]; then
        echo "Error in bedtools intersect for $base_name!" >> $LOGFILE
        exit 1
    fi

    # Counting remaining reads
    echo "Counting remaining reads for sample $base_name..." >> $LOGFILE
    remaining_reads=$(samtools view -c -F 260 ${base_name}.bam)
    echo "Remaining reads for sample $base_name: $remaining_reads" >> $LOGFILE

    # Store the remaining reads to determine the lowest value later for splitting
    remaining_reads_array+=($remaining_reads)
    
    # Indexing bam file
    samtools index -@ $THREADS -b $base_name.bam

    # Creating bedgraph file
    echo "Creating bedgraph for sample $base_name..." >> $LOGFILE
    bedtools genomecov -bg -pc -ibam ${base_name}.bam > ${base_name}.bedgraph
    
    source activate $DEEPTOOLS
    
    # Creating bigwig files
    echo "Creating bigwig for sample $base_name..." >> $LOGFILE
    bamCoverage -b ${base_name}.bam -o ${base_name}.bw -bs 10 -p $THREADS --normalizeUsing RPGC --effectiveGenomeSize 2701495711 -e >> $LOGFILE 2>&1

    # Removing intermediate files
    rm ${base_name}*.fastq
    rm ${base_name}.bowtie2.*

done

source activate $CUTNRUN

# Find the lowest remaining reads count from all replicates
min_reads=$(printf "%s\n" "${remaining_reads_array[@]}" | sort -n | head -n 1)
echo "Lowest number of remaining reads: $min_reads" >> $LOGFILE

# Loop through BAM files and split into aggregates
for sample in *_1.bam; do
    for i in $(seq 1 $REPLICATES); do
        for aggregate in Aggregate1 Aggregate2 Aggregate3; do
            echo "Shuffling and splitting sample for replicate $i into $aggregate..." >> $LOGFILE
            # Dynamically set the split size to the lowest remaining reads count
            samtools view ${sample%_1.bam}_${i}.bam | shuf - | split -d -l $min_reads - ${aggregate}/${sample%_1.bam}_${i}
        done
    done
done

# Build aggregates and call peaks with MACS2
for aggregate in Aggregate1 Aggregate2 Aggregate3; do
    source activate $CUTNRUN
    echo "Building aggregate for $aggregate..." >> $LOGFILE
    cat header.sam ${aggregate}/*_100 | samtools view -@$THREADS -bS | samtools sort -@$THREADS > ${aggregate}.bam
    
    # Perform peak calling for the first aggregate (without replicates)
    macs2 callpeak -t ${aggregate}.bam -c $IgG -f BAMPE --keep-dup all -n ${aggregate}/${aggregate}_1 $MACS2_THRESH >> $LOGFILE 2>&1
    peak_count_1=$(wc -l ${aggregate}/${aggregate}_1_peaks.narrowPeak | awk '{print $1}')
    echo "$peak_count_1" >> ${aggregate}.txt
    
    cp ${aggregate}.bam temp2.bam
    
    # Initialize a variable to track the last replicate added
    last_replicate=""
    
    # Loop through replicates for peak calling
    for i in $(seq 2 $REPLICATES | shuf); do
        echo "Processing replicate $i for aggregate $aggregate..." >> $LOGFILE
        cat header.sam ${aggregate}/*_${i}00 | samtools view -@$THREADS -bS | samtools sort -@$THREADS > temp_sort.bam
        samtools merge ${aggregate}.bam temp2.bam temp_sort.bam -@$THREADS -f
        cp ${aggregate}.bam temp2.bam
        
        # Perform peak calling after adding each replicate
        macs2 callpeak -t ${aggregate}.bam -c $IgG -f BAMPE --keep-dup all -n ${aggregate}/${aggregate}_${i} $MACS2_THRESH >> $LOGFILE 2>&1
        
        # Count the peaks after this replicate is added
        peak_count=$(wc -l ${aggregate}/${aggregate}_${i}_peaks.narrowPeak | awk '{print $1}')
        echo "$peak_count" >> ${aggregate}.txt
        
        # Update last_replicate to the most recent replicate
        last_replicate=$i
    done
    
    # Final output: Name the file with the last replicate added
    final_output="${aggregate}/${aggregate}_final_peaks.narrowPeak"
    cp ${aggregate}/${aggregate}_${last_replicate}_peaks.narrowPeak $final_output
    echo "Final output for $aggregate saved as: $final_output"
    
    samtools index -@ $THREADS -b ${aggregate}.bam
    
    # Make bigwig of each aggregate for visualization
    source activate $DEEPTOOLS
    echo "Creating bigwig for aggregate $aggregate..." >> $LOGFILE
    bamCoverage -b ${aggregate}.bam -o ${aggregate}.bw -bs 10 -p $THREADS --normalizeUsing RPGC --effectiveGenomeSize 2701495711 -e 
    
done

# Make average bigwig of three aggregates for visualization.
source activate $DEEPTOOLS
echo "Creating ICEBERG bigwig" >> $LOGFILE 
bigwigAverage -p $THREADS -b Aggregate1.bw Aggregate2.bw Aggregate3.bw -o ICEBERG.bw >> $LOGFILE 2>&1

source activate $CUTNRUN

# Define ICEBERG peaks
echo "Defining ICEBERG peaks..." >> $LOGFILE
bedtools intersect -wa -a Aggregate1/Aggregate1_final_peaks.narrowPeak -b Aggregate2/Aggregate2_final_peaks.narrowPeak > del1.bed
bedtools intersect -wa -a Aggregate2/Aggregate2_final_peaks.narrowPeak -b Aggregate1/Aggregate1_final_peaks.narrowPeak > del2.bed
bedtools intersect -wa -a Aggregate1/Aggregate1_final_peaks.narrowPeak -b Aggregate3/Aggregate3_final_peaks.narrowPeak > del3.bed
bedtools intersect -wa -a Aggregate3/Aggregate3_final_peaks.narrowPeak -b Aggregate1/Aggregate1_final_peaks.narrowPeak > del4.bed
bedtools intersect -wa -a Aggregate2/Aggregate2_final_peaks.narrowPeak -b Aggregate3/Aggregate3_final_peaks.narrowPeak > del5.bed
bedtools intersect -wa -a Aggregate3/Aggregate3_final_peaks.narrowPeak -b Aggregate2/Aggregate2_final_peaks.narrowPeak > del6.bed
cat del*.bed > del7.bed
bedtools sort -i del7.bed > del8.bed
bedtools merge -i del8.bed > ICEBERG_Peaks.bed

# Log number of final peaks
wc -l ICEBERG_Peaks.bed >> $LOGFILE

# Clean up temporary files
rm del*.bed
rm temp*.bam

# Make heatmaps of all replicates and ICEBERG bw within ICEBERG peaks
source activate $DEEPTOOLS
echo "Creating heatmap for ICEBERG peaks" >> $LOGFILE
computeMatrix reference-point -R ICEBERG_Peaks.bed -S ICEBERG.bw *_*.bw -o matrix -b 2000 -a 2000 -p $THREADS --referencePoint center -bs 50 --missingDataAsZero >> $LOGFILE 2>&1
plotHeatmap -m matrix -o Heatmap_ICEBERG.pdf --colorMap Greys >> $LOGFILE 2>&1

# Process aggregate.txt files to create a plot-compatible format
echo "Replicate PeakCount" > Aggregate1_plot_data.txt
awk '{print NR, $1}' Aggregate1.txt >> Aggregate1_plot_data.txt

echo "Replicate PeakCount" > Aggregate2_plot_data.txt
awk '{print NR, $1}' Aggregate2.txt >> Aggregate2_plot_data.txt

echo "Replicate PeakCount" > Aggregate3_plot_data.txt
awk '{print NR, $1}' Aggregate3.txt >> Aggregate3_plot_data.txt

# Now, create the plot using gnuplot
echo "
set terminal pdf size 6,4
set output 'aggregate_peaks.pdf'
set title 'Peak Counts Across Replicates'
set xlabel 'Replicates'
set ylabel 'Peak Count'
set grid
set key left top

plot 'Aggregate1_plot_data.txt' using 1:2 with linespoints title 'Aggregate 1', \
     'Aggregate2_plot_data.txt' using 1:2 with linespoints title 'Aggregate 2', \
     'Aggregate3_plot_data.txt' using 1:2 with linespoints title 'Aggregate 3'
" > plot_script.gp

# Run the plot script with gnuplot
gnuplot plot_script.gp

# Clean up the temporary files
rm plot_script.gp Aggregate1_plot_data.txt Aggregate2_plot_data.txt Aggregate3_plot_data.txt

# Final message
echo "Pipeline execution completed: $(date)" >> $LOGFILE
echo "ICEBERG pipeline completed successfully :D :D :D"
