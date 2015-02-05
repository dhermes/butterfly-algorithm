help:
	@echo 'Makefile for Butterfly Slides                         '
	@echo '                                                      '
	@echo 'Usage:                                                '
	@echo '   make pdf        create Beamer PDF from LaTeX       '
	@echo '   make slides     render slides from IPython Notebook'

pdf:
	pdflatex butterfly_slides.tex

slides:
	ipython nbconvert butterfly.ipynb --to slides

.PHONY: help pdf slides
