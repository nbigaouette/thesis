Wail Gueaieb and Stephen M. Carr
Last updated: 2005-04-22

This "ready-made thesis template" is meant to help the University of Ottawa's
graduate students to write their theses in LaTeX. It is a modified version of
the package made by Stephen M. Carr for the University of Waterloo (Thank you
Stephen). It is NOT an official release of the University of Ottawa.
Nevertheless, it is made to satisfy the thesis' technical standards set by the
University of Ottawa at:
http://www.grad.uottawa.ca/regulations/thesis_research/manual/appendixb.html
http://www.grad.uottawa.ca/regulations/thesis_research/manual/appendixc.html

The template was originally prepared for SITE students. Students in other 
disciplines may also use it, however they should double check with their 
academic units. Remember that it is the student's responsibility to make sure 
that sources are referenced properly. The Faculty of Graduate and Postdoctoral 
studies at the University of Ottawa mendates that quotation marks are used in 
direct quotations.

This package comes with NO WARRANTY and NO SUPPORT! You may use it at YOUR OWN
RISK!!!

Package Content
---------------
bibliography:
   keylatex.bib:
      This is a bibliography file that contains the references to be inserted
in the Bibliography part of the thesis.

private:
   thesis-front-pages.tex:
   thesis-headings.tex:
   thesis-preamble.tex:
      These 3 files contain the commands necessary to conform to the thesis'
technical standards of the University of Ottawa. Do not edit them unless you
know what you are doing.

00readme.txt:
   This file.

abstract.tex:
   Contains the abstract of the thesis.

acknowledgement.tex:
   Contains the acknowledgement of the thesis.

appendix1-help.tex:
   Contains the content of Appendix 1 of the thesis.

beam.eps:
   A figure.

chapterN-XXXX.tex:
   Chapter N of the thesis.

glossary.tex:
   Thesis' glossary.

mastersource.tex:
   A file used in chapter2-creating.tex.

uottawa-thesis.tex:
   The main file of the thesis.


Command Sequence Used to Process These Files
--------------------------------------------
At the Unix prompt (tetex distribution) or DOS / Command Shell (MikTeX
distribution under Windows 95/98/NT/2000), process these source files
with the following command sequence:

latex uottawa-thesis
bibtex uottawa-thesis
makeindex uottawa-thesis
latex uottawa-thesis
latex uottawa-thesis

This much creates the "dvi" file, uottawa-thesis.dvi, which can be previewed.
To create a printable Postscript file, post-process the "dvi" file with

dvips uottawa-thesis.dvi > uottawa-thesis.ps

To create an Adobe PDF file for distribution electronically,

dvips -Ppdf uottawa-thesis.dvi > uottawa-thesis.ps
This uses the Type 1 Postcript versions of the Computer Modern fonts.
Then convert the PS file into PDF with Adobe Acrobat 4.0 (Writer), or with
the freeware Ghostview 4.3 program (using File/Convert, then select pdfwrite
and 600 dpi resolution).

Alternately, use the "dvipdfm" program (distributed with MikTeX) to postprocess
the "dvi" file directly into PDF:

dvipdfm uottawa-thesis.dvi