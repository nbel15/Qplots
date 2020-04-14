# Qplots
QPlots is a user-friendly interface built under R language using Shiny package, it allows to rapidly explore microbiome data

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
