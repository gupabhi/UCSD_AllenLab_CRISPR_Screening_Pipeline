Sweave("Nitrogen__countsummary.Rnw");
library(tools);

texi2dvi("Nitrogen__countsummary.tex",pdf=TRUE);

