all: stripe82.pdf

%.pdf: %.tex
	pdflatex $<
	bash -c " ( grep 'undefined citations' $*.log && bibtex $*.aux; pdflatex $*.tex ) || echo noRerun "
	bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "

