help:
	@echo 'Makefile for Butterfly Slides                                 '
	@echo '                                                              '
	@echo 'Usage:                                                        '
	@echo '   make pdf                create Beamer PDF from LaTeX       '
	@echo '   make slides             render slides from IPython Notebook'
	@echo '   make serve-slides       serve slides from IPython Notebook '
	@echo '   make _fortran_utils.so  build Fortran speedup module       '

pdf:
	pdflatex butterfly_slides.tex

slides:
	MAKE_BUTTERFLY_SLIDES="True" ipython nbconvert butterfly.ipynb --to slides

serve-slides:
	MAKE_BUTTERFLY_SLIDES="True" ipython nbconvert butterfly.ipynb --to slides --post serve

_fortran_utils.so:
	f2py -c -m _fortran_utils speedup.f90

.PHONY: help pdf slides serve-slides
