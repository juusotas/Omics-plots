# Omics-plots
Dump for different functions to plot gene set analysis data.

This repository holds all the ploting functions I use in my gene set analysis pipeline.
These functions are used through the plot_gsea_results functions, which loops through a list of gene set analysis results. These analysis results are then used to plot different types of plots, which some have more and others less utility depending on the gsea results.
For example, the violin_plotly functions plots a interactive html violin plot, which allows you to study individual genes and their fold changes in n-pathways, this way is especially great if only a few pathways have change, but becomes less informatic if multiple have changed. A more usefull plot in the case of multiple significantly altered pathways is the gsea volcano plot, which plots NES on the X-axis and Pathways on the Y-axis, the unfortunate thing with this plot is that you lose the gene level information, which could also be very informative.
