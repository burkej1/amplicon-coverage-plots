# Takes an output file as first argument
BASENAME=${1}

SUFFIXARRAY=( "_sample-normalised_coverage.html" "_amplicon-normalised_coverage.html" "_all.html" )

# Creating these separately due to regex matching too much
DEPTHPLOTS=$(ls coverage_plots/* | grep _coverage.html | grep -ve sample-normalised_coverage -ve amplicon-normalised_coverage)
DEPTHOUTPUT="${BASENAME}_coverage.html"
echo ${DEPTHOUTPUT}
# Top bit of the html
echo "<!DOCTYPE html>" > ${DEPTHOUTPUT}
echo "<html>" >> ${DEPTHOUTPUT}
echo "<body>" >> ${DEPTHOUTPUT}
# Add iframes for each plot assuming plots are in a folder called "coverage_plots"
for HTML in ${DEPTHPLOTS[@]}; do
    echo "<iframe src='${HTML}' style='border:none;' width='100%' height='1000'>" >> ${DEPTHOUTPUT}
    echo "</iframe>" >> ${DEPTHOUTPUT}
done
# Add the end
echo "</body>" >> ${DEPTHOUTPUT}
echo "</html>" >> ${DEPTHOUTPUT}

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
