all: Sim.pdf
Sim.pdf: Sim.tex \
         Atmosphere.tex \
         LGS.tex \
         SHWFS.tex \
         DM.tex \
         Control.tex \
         Derivatives.tex \
         DF.tex \
         Fresnel.tex \
         Biblio.tex
	pdflatex --file-line-error-style Sim.tex
	makeindex Sim.idx
	pdflatex --file-line-error-style Sim.tex

clean:
	rm -f *.aux *.out *.log *.ilg *.idx *.ind *.toc