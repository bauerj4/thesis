#!/bin/bash

files=("myPreamble" "Abstract" "Introduction" "Background" "Paper_I" "Implementation" "Paper_II" "Paper_III" "Conclusion" "Appendices") 
#files[0]="myPreamble.tex" 
#files[1]="Abstract.tex" 
#files[2]="Background.tex" 
#files[3]="Paper_I.aux" 
#files[4]="Conclusion.tex" 
pdflatex Thesis.tex
for auxfile in "${files[@]}"
do
    echo "Opening file"
    bibtex `basename $auxfile .aux`
done
pdflatex Thesis.tex
pdflatex Thesis.tex
