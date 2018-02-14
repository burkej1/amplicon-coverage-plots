# Takes an output file as first argument
BASENAME=${1}

SUFFIXARRAY=( "_sample_normalised.html" "_amplicon_normalised.html" "_coverage.html" "_all.html" )

for SUFFIX in ${SUFFIXARRAY[@]}; do
    echo ${BASENAME}${SUFFIX}
    # Top bit of the html
    echo "<!DOCTYPE html>" > ${BASENAME}${SUFFIX}
    echo "<html>" >> ${BASENAME}${SUFFIX}
    echo "<body>" >> ${BASENAME}${SUFFIX}
    # Add iframes for each plot assuming plots are in a folder called "coverage_plots"
    for plot in coverage_plots/*${SUFFIX}; do
        echo "<iframe src='${plot}' style='border:none;' width='100%' height='1000'>" >> ${BASENAME}${SUFFIX}
        echo "</iframe>" >> ${BASENAME}${SUFFIX}
    done
    # Add the end
    echo "</body>" >> ${BASENAME}${SUFFIX}
    echo "</html>" >> ${BASENAME}${SUFFIX}
done
