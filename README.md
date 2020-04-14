# Qplots
QPlots is a user-friendly interface built under R language using Shiny package, it allows to rapidly explore microbiome data

For running Qplot app:
1. Install required packages using "install_packages.R" file
2. Open "server.R" file and click on "Run App" button in top right of the script window (or Ctrl + Shift + Enter)

Description of the content
Folders:
- "analysis_functions": contains the main functions of the downstream analysis (relative abundance, Alpha diversity, ....)
- "server_functions": contains the server logic required for the analysis sections
- "www": Contains the guide sections in HTML format and the example files

Files:
- server.R: contains the server function definition, it is require the files saved in "server_functions" folder
- ui.R: contains the user interface definition
- tools_functions: contains six extra functions needed during the analysis
- fun.R: contains additional functions needed for the data preparation 
