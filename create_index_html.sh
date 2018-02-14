# Takes an output file as first argument
rm ${1}

# Top bit of the html
echo "<!DOCTYPE html>" >> ${1}
echo "<html>" >> ${1}
echo "<body>" >> ${1}

# Add iframes for each plot assuming plots are in a folder called "coverage_plots"
for plot in coverage_plots/*.html; do
    echo "<iframe src='${plot}' style='border:none;' width='100%' height='1000'>" >> ${1}
    echo "</iframe>" >> ${1}
done

# Add the end
echo "</body>" >> ${1}
echo "</html>" >> ${1}
