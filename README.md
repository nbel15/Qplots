# Qplots
Qplots is a Shiny-based web application for interactive analysis and visualization of microbial sequencing data. Qplots constitute an analytical pipeline for the calculation of the bacterial composition and diversity based on Operational Taxonomic Units (OTUs) tables. Currently, it provides a set of well-documented choices of downstream analysis including taxonomic composition, alpha and beta diversity analysis, statistical comparisons between subsets of samples, as well as a variety of visualization methods including tables, boxplots, heatmaps, bar charts, and different ordination plots. Qplots was implemented entirely in R language using the Shiny framework. It can be easily used locally in any system with R installed, including Windows, macOS, and most Linux distributions, or remotely through a web server without bioinformatics expertise. It can also be used as a framework for advanced users who can modify and expand the tool. With Qplots, we aim to provide the scientific community with a platform that simplifies the analysis and enables publication-quality visualization of microbiome data.
![image](https://user-images.githubusercontent.com/19431517/169152457-f858b15d-1c8d-4354-bfbe-9c81d599da06.png)


<h3>For running Qplot app: </h3>
<ol type="1">
  <li> Install required packages using "install_packages.R" file</li>
  <li> Open "server.R" file and click on "Run App" button in top right of the script window (or Ctrl + Shift + Enter)</li>
</ol>

<h3>Description of the content</h3>

<h4>Folders:</h4>
<ul>
  <li> "analysis_functions": contains the main functions of the downstream analysis (relative abundance, Alpha diversity, ....)</li>
  <li> "server_functions": contains the server logic required for the analysis sections</li>
  <li> "www": Contains the guide sections in HTML format and the example files</li>
</ul>

<h4>Files:</h4>
<ul>
  <li> server.R: contains the server function definition, it is require the files saved in "server_functions" folder </li>
  <li> ui.R: contains the user interface definition</li>
  <li> tools_functions: contains six extra functions needed during the analysis</li>
  <li> fun.R: contains additional functions needed for the data preparation </li>
</ul>
