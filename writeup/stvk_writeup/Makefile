# type "make" command in Unix to create asme2ej.pdf file
all:
	latex writeup
	bibtex writeup
	latex writeup
	latex writeup
	dvipdfm writeup

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
