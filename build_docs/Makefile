# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
APIDOCDIR     = $(SOURCEDIR)/auto
BUILDDIR      = _build
AUTOSUMMARYDIR= _autosummary
PUBDIR        = ../docs
MODULESDIR    = ../src/bgc_md2

# Internal variables.
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(SPHINXOPTS)

# Put it first so that "make" without argument is like "make help".
.PHONY: help 
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/*
	rm -rf $(APIDOCDIR)/*
	rm -rf $(SOURCEDIR)/$(AUTOSUMMARYDIR)/*

.PHONY: api-doc
api-doc:
	sphinx-apidoc -o $(APIDOCDIR)/ --ext-autodoc $(MODULESDIR)

.PHONY: html
html: api-doc
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

.PHONY: publish 
publish:
	@rm -rf $(PUBDIR)/* 
	@cp -r $(BUILDDIR)/html/* $(PUBDIR)
	@touch $(PUBDIR)/.nojekyll #tell github to serve plain html


.PHONY: help Makefile


# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
