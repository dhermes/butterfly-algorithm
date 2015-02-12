help:
	@echo 'Makefile for Butterfly Slides                                   '
	@echo '                                                                '
	@echo 'Usage:                                                          '
	@echo '   make slides               render slides from IPython Notebook'
	@echo '   make serve-slides         serve slides from IPython Notebook '
	@echo '   make serve-static-slides  serve slides from HTML             '
	@echo '   make generated-module     make module w/methods from Notebook'

slides:
	MAKE_BUTTERFLY_SLIDES="True" ipython nbconvert butterfly.ipynb --to slides

serve-slides:
	ipython nbconvert butterfly.ipynb --to slides --post serve

serve-static-slides:
	cd slides-static/ && python -m SimpleHTTPServer

generated-module:
	python only_functions_from_ast.py

.PHONY: help pdf slides serve-slides serve-static-slides generated-module
