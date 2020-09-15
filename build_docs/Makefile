# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      = _build
AUTOSUMMARYDIR= _autosummary
PUBDIR        = ../docs

# Internal variables.
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(SPHINXOPTS)

# Put it first so that "make" without argument is like "make help".
.PHONY: help 
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/*
	rm -rf $(SOURCEDIR)/$(AUTOSUMMARYDIR)/*

.PHONY: html
html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."



.PHONY: publish 
publish:
	@rm -rf $(PUBDIR)/* 
	@cp -r $(BUILDDIR)/html/* $(PUBDIR)
	@touch $(PUBDIR)/.nojekyll #tell github to serve plain html

