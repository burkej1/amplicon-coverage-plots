INTERVALS="/Users/Jared/Documents/Code/amplicon-coverage-plots/testing_files/BRA-STRAP_621717_100.final.roverfile_g37.bed"
OUTPUTFILE="amplicon_metrics.tsv"

# Remove any leftover files from previous runs
rm ${OUTPUTFILE}

for bam in testing_files/*.bam; do
    SAMPLE=${bam%.clipped.sort.hq.bam}  # Get the sample string
    SAMPLE=${SAMPLE##*/}  # Remove the path
    bedtools coverage -f 5E-1 -a ${INTERVALS} -b ${bam} | sed "s/$/	${SAMPLE}/g" >> ${OUTPUTFILE} # Calculate amplicon coverage metrics
done

