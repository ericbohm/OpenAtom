all: clean fullpdf fullhtml flowpdf flowhtml

fullpdf: OpenAtomFullDoc.tex
	rm -r OpenAtomFullDoc | true
	pdflatex OpenAtomFullDoc.tex

fullhtml: OpenAtomFullDoc.tex
	rm -f *.pdf
	rm -r OpenAtomFlow | true
	latex2html -local_icons OpenAtomFullDoc.tex

flowpdf: OpenAtomFlow.tex
	rm -f *.pdf
	pdflatex OpenAtomFlow.tex 
	pdflatex OpenAtomFlow.tex 
	bibtex OpenAtomFlow
	pdflatex OpenAtomFlow.tex 

flowhtml: OpenAtomFlow.tex
	latex2html -local_icons OpenAtomFlow.tex
	latex2html -local_icons OpenAtomFlow.tex
	bibtex OpenAtomFlow
	latex2html -local_icons OpenAtomFlow.tex

clean:
	rm -f *.ilg *.aux *.log *.dvi *.idx *.toc *.lof *.lot
	rm -f *.blg *.bbl *~

