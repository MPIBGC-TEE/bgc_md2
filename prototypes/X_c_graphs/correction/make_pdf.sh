#!/bin/bash
f=manusscript.tex
#rm *.blg
#rm *.log 
#rm *.aux
pdflatex $f
bibtex
pdflatex $f
pdflatex $f
