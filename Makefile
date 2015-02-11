help:
	@echo 'Makefile for Butterfly Slides                                   '
	@echo '                                                                '
	@echo 'Usage:                                                          '
	@echo '   make pdf                  create Beamer PDF from LaTeX       '
	@echo '   make slides               render slides from IPython Notebook'
	@echo '   make serve-slides         serve slides from IPython Notebook '
	@echo '   make serve-static-slides  serve slides from HTML             '

pdf:
	pdflatex butterfly_slides.tex

slides:
	MAKE_BUTTERFLY_SLIDES="True" ipython nbconvert butterfly.ipynb --to slides

serve-slides:
	ipython nbconvert butterfly.ipynb --to slides --post serve

serve-static-slides:
	python -m SimpleHTTPServer

.PHONY: help pdf slides serve-slides
