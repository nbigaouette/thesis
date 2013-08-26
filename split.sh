#!/bin/bash

# Color pages
color_pages=(7 10 13 14 15 18 21 31 35 36 39 43 96 97 98 99 100 101 112 113)
# Shift the pages since title page to table of content have roman numbers
nb_roman_pages=13
for ((i = 0 ; i < ${#color_pages[*]} ; i++)); do
    color_pages_real[${i}]=$((${color_pages[${i}]} + ${nb_roman_pages}))
done

filename_color="NBigaouette_Thesis_Colour.pdf"
filename_bw="NBigaouette_Thesis_BW.pdf"
filename_basis="thesis.pdf"

echo "Colour pages: ${color_pages[*]} (${color_pages_real[*]})"

rm -f ${filename_color} ${filename_bw}

stapler sel ${filename_basis} ${color_pages_real[*]} ${filename_color}
stapler del ${filename_basis} ${color_pages_real[*]} ${filename_bw}

# # Convert B/W pdf to real greyscale
# # See http://superuser.com/questions/104656/convert-a-pdf-to-greyscale-on-the-command-line-in-floss
# gs \
#  -sOutputFile=output.pdf \
#  -sDEVICE=pdfwrite \
#  -sColorConversionStrategy=Gray \
#  -dProcessColorModel=/DeviceGray \
#  -dCompatibilityLevel=1.4 \
#  -dNOPAUSE \
#  -dBATCH \
#  -dAutoRotatePages=/None \
#  ${filename_bw}

echo "Done!"
