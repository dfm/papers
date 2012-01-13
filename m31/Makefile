.SUFFIXES:
.SUFFIXES: .tex .pdf .ps

all: m31.pdf

.tex.pdf:
	pdflatex $*.tex
	- bash -c " ( grep 'undefined citations' $*.log && bibtex $*.aux; pdflatex $*.tex ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $*.tex ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $*.tex ) || echo noRerun "

.tex.ps:
	latex $*.tex
	- bash -c " ( grep Rerun $*.log && latex $*.tex ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && latex $*.tex ) || echo noRerun "
	- bash -c " dvips $*.dvi "

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.dvi *.pdf
