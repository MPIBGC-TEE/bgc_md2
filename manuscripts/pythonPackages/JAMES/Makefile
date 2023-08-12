.PHONY: clean all spellcheck
all: main.tex TEE-clean.bib 
	pdflatex main.tex
	#bibexport -o minimal.bib main.aux
	bibtex main.aux
	pdflatex main.tex
	pdflatex main.tex

clean:
	rm *.aux *.bbl *.blg
spellcheck:
	aspell -c -t content.tex -d en
