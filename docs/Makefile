# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    =
SOURCEDIR      = source
BUILDDIR      = _build

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don't have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS)
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  latexpdfja to make LaTeX files and run them through platex/dvipdfmx"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  xml        to make Docutils-native XML files"
	@echo "  pseudoxml  to make pseudoxml-XML files for display purposes"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"

clean:
	rm -rf $(BUILDDIR)/*

html:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/html $(ALLSPHINXOPTS) -b html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/dirhtml $(ALLSPHINXOPTS) -b dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/singlehtml $(ALLSPHINXOPTS) -b singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/pickle $(ALLSPHINXOPTS) -b pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/json $(ALLSPHINXOPTS) -b json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/htmlhelp $(ALLSPHINXOPTS) -b htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/qthelp $(ALLSPHINXOPTS) -b qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/mbuild.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/mbuild.qhc"

devhelp:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/devhelp $(ALLSPHINXOPTS) -b devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/mbuild"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/mbuild"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/epub $(ALLSPHINXOPTS) -b epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/latex $(ALLSPHINXOPTS) -b latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/latex $(ALLSPHINXOPTS) -b latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

latexpdfja:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/latex $(ALLSPHINXOPTS) -b latex
	@echo "Running LaTeX files through platex and dvipdfmx..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf-ja
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/text $(ALLSPHINXOPTS) -b text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/man $(ALLSPHINXOPTS) -b man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/texinfo $(ALLSPHINXOPTS) -b texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/texinfo -b texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/locale $(I18NSPHINXOPTS)  -b gettext
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/changes $(ALLSPHINXOPTS) -b changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/linkcheck $(ALLSPHINXOPTS) -b linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/doctest $(ALLSPHINXOPTS) -b doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

xml:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/xml $(ALLSPHINXOPTS) -b xml
	@echo
	@echo "Build finished. The XML files are in $(BUILDDIR)/xml."

pseudoxml:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)/pseudoxml $(ALLSPHINXOPTS) -b pseudoxml
	@echo
	@echo "Build finished. The pseudo-XML files are in $(BUILDDIR)/pseudoxml."
