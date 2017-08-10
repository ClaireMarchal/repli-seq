your_genome_chrom_sizes="" # Paste between the quotes the path to your genome index file, including the file name.
path_to_your_genome="" # Paste between the quotes the path to your bowtie genome indexes including the prefix of the files.
size=50000 # Change the windows size here if needed (size is in bp).

# Generating 50kb windows positions along the genome. You can change the -w value to change the windows size, and the -s value, to change the steps (to make overlapping windows for example)
sort -k1,1 -k2,2n $your_genome_chrom_sizes > your_genome_sorted.chrom.sizes
bedtools makewindows -w $size -s $size -g your_genome_sorted.chrom.sizes > your_genome_windows.bed

# Mapping fastq(.gz) and making bed with values on genomic windows
for file in *.fastq*; do
    # Mapping
    bowtie2 -x $path_to_your_genome --no-mixed --no-discordant --reorder -U $file -S ${file%.fastq*}.sam 2>> ${file%.fastq*}_mapping_log.txt
    # sam to bam conversion
    samtools view -bSq 20 ${file%.fastq*}.sam > ${file%.fastq*}.bam
    # duplicate reads suppression
    samtools sort -o ${file%.fastq*}_srt.bam ${file%.fastq*}.bam
    samtools rmdup -S ${file%.fastq*}_srt.bam ${file%.fastq*}_rmdup.bam
    # bam to bed conversion
    bamToBed -i ${file%.fastq*}_rmdup.bam | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S 5G > ${file%.fastq*}.bed
    # bed line number calculation
    x=`wc -l ${file%.fastq*}.bed | cut -d' ' -f 1`
    # generate coverage on genomic windows
    bedtools intersect -sorted -c -b ${file%.fastq*}.bed -a your_genome_windows.bed | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > ${file%.fastq*}.bg
done

# Calculating RT
n=0
for file in *_E_.bg; do
    paste $file ${file%E_.bg}L_.bg | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > ${file%E_.bg}T_.bg
    n=$((n+1));
done

# Merging RT files
if [ "$n" -gt "1" ]; then
    bedtools unionbedg -filler "NA" -i *T_.bg > merge_RT.txt
else
    cp *T_.bg merge_RT.txt
fi

