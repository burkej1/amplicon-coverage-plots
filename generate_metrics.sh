INTERVALS="/Users/Jared/Documents/Code/amplicon-coverage-plots/testing_files/BRA-STRAP_621717_100.final.roverfile_g37.bed"
OUTPUTFILE="amplicon_metrics.tsv"

# Remove any leftover files from previous runs
rm ${OUTPUTFILE}
rm .temp-coverage
rm .temp-samples

for bam in testing_files/*.bam; do
    bedtools coverage -a ${INTERVALS} -b ${bam} > .temp-coverage  # Calculate amplicon coverage metrics
    SAMPLE=${bam%.clipped.sort.hq.bam}  # Get the sample string
    SAMPLE=${SAMPLE##*/}  # Remove the path
    LINENUMBER=$(wc -l <${INTERVALS})  # Get a line number (for generating the sample column)
    # Loop to create a column containing the sample ID repeated
    while [ ${LINENUMBER} -ne 0 ]; do
        echo ${SAMPLE} >> .temp-samples
        LINENUMBER=$[${LINENUMBER} - 1]
    done
    paste .temp-coverage .temp-samples >> ${OUTPUTFILE}  # Paste the sample column
    rm .temp-samples  # Emptying file to prepare for next loop
done

rm .temp-coverage
rm .temp-samples
