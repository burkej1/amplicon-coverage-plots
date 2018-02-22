# Takes an output file as first argument
BASENAME=${1}
DIRECTORY=${2%/}

DEPTHPLOTS=$(ls ${DIRECTORY}/*)
DEPTHOUTPUT="${BASENAME}_coverage.html"

echo ${DEPTHOUTPUT}
# Top bit of the html
echo "<!DOCTYPE html>" > ${DEPTHOUTPUT}
echo "<html>" >> ${DEPTHOUTPUT}
echo "<body>" >> ${DEPTHOUTPUT}
# Add iframes for each plot
for HTML in ${DEPTHPLOTS[@]}; do
    echo "<iframe src='${HTML}' style='border:none;' width='100%' height='1000'>" >> ${DEPTHOUTPUT}
    echo "</iframe>" >> ${DEPTHOUTPUT}
done
# Add the end
echo "</body>" >> ${DEPTHOUTPUT}
echo "</html>" >> ${DEPTHOUTPUT}


# # Ignore this
# SUFFIXARRAY=( "_sample-normalised_coverage.html" "_amplicon-normalised_coverage.html" "_all.html" )

# for SUFFIX in ${SUFFIXARRAY[@]}; do
#     echo ${BASENAME}${SUFFIX}
#     # Top bit of the html
#     echo "<!DOCTYPE html>" > ${BASENAME}${SUFFIX}
#     echo "<html>" >> ${BASENAME}${SUFFIX}
#     echo "<body>" >> ${BASENAME}${SUFFIX}
#     # Add iframes for each plot
#     for plot in coverage_plots/*${SUFFIX}; do
#         echo "<iframe src='${plot}' style='border:none;' width='100%' height='1000'>" >> ${BASENAME}${SUFFIX}
#         echo "</iframe>" >> ${BASENAME}${SUFFIX}
#     done
#     # Add the end
#     echo "</body>" >> ${BASENAME}${SUFFIX}
#     echo "</html>" >> ${BASENAME}${SUFFIX}
# done
