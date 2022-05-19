library(shinydashboard)
library(shiny)
library(varhandle)
library(gplots)
library(ggplot2)
library(shinyBS)
library(shinyjs)
library(shinythemes)
library(plotly)
library(reshape2)
library(DT)
library(magrittr)
library(dplyr)
library(RColorBrewer)
library(networkD3)
library(igraph)
library(tidyr)
library(ade4)
library(GUniFrac)
library(phangorn)
library(cluster)
library(fpc)
library(assertr)
library(units)
library(ape)
library(picante)

library(BiocManager)

library(adespatial)
library(vegan)

library(ggtree)
library(phyloseq)
library(viridis)
library(viridis)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)


library(QsRutils)

library(shinydashboardPlus)

library(ggpubr)
library(pairwiseAdonis)
library(tidyr)
library(stringr)
library(shinyalert)

library(mattsUtils)
library(ggvenn)
library(microbiome)

dashboardPage(
  dashboardHeader(
    title = "Amplicon analysis"
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Guide", tabName = "guide", icon = icon("book-reader")),
      menuItem("Import Data", tabName = "import_data", icon = icon("file-import")),
      menuItem("Beta diversity", tabName = "beta_div_tab", icon = icon("ruler-combined"),
               menuSubItem("Graphs", tabName = "beta_div_graph")
               #menuSubItem("Statistical Analysis", tabName = "beta_div_stat")
               ),
      menuItem("Alpha diversity", tabName = "alpha_tab", icon = icon("project-diagram")),
      menuItem("Relative Abundance", tabName = "rb_tab", icon = icon("poll")),
      menuItem("Core microbiome", tabName = "core_tab", icon = icon("linode"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    useShinyalert(),
    
    tabItems(
      tabItem(tabName = "guide",
              fluidRow(
                        tabBox(
                          width = 12, 
                          tabPanel("Qplots Workflow", uiOutput("qplot_workflow")),
                          tabPanel("FastQ analysis", uiOutput("script_guide"))
                        )) 
      ),
      tabItem(tabName = "import_data", 
              fluidRow(
                box(
                  width = 12,
                  title = "Load Data", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  column(width = 3, 
                         fileInput('otumat_input', 'Choose OTUs file',
                                   accept=c('text', 
                                            'text/comma-separated-values,text/plain'))
                         ),
                  column(width = 3, 
                         fileInput('taxamat_input', 'Choose Taxa file',
                                   accept=c('text/tsv', 
                                            'text/comma-separated-values,text/plain',
                                            '.tsv'))),
                  column(width = 3, 
                         fileInput('metamat_input', 'Choose mapping file',
                                   accept=c('text', 
                                            'text/comma-separated-values,text/plain'))
                         ),
                  column(width = 3, 
                         fileInput('tree_file_input', 'Choose tree file',
                                   accept=c(".nwk",'.NWK','.tre', '.tree'))
                         
                         ),
                  actionButton("read_input", "Read Input")
                  )
                ),
              fluidRow(
                tabBox(width = 12, id = "tabview",
                       title = "Visualize Input",
                       tabPanel("OTUs table",
                                div(style = 'overflow-x: scroll', DT::dataTableOutput("otu_tab_view",width = "100%"))),
                       tabPanel("Taxonomy table", 
                                div(style = 'overflow-x: scroll', DT::dataTableOutput("taxa_tab_view",width = "100%"))),
                       tabPanel("Metatable", 
                                div(style = 'overflow-x: scroll', DT::dataTableOutput("meta_tab_view",width = "100%")))
                       )
                )
              ),
      tabItem(tabName = "beta_div_graph",
              box(
                width = 12, title = "Beta diversity parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                fluidRow(
                  column(width = 3, 
                         selectInput("beta.div.factor", "Factor", choices = "choose")
                         ),
                  column(width = 3,
                         selectInput("dist.method", "Distance metric",
                                     choices = c("GUnifrac" = "d_0.5", "Weighted" = "d_1", "Unweighted" = "d_UW", 
                                                 "VAW Unifract" = "d_VAW", "Bray Curtis"="bray"), selected = "d_0.5")
                         ),
                  column(width = 3,
                         selectInput("ordination.method", "Ordination", 
                                     choices = c("MDS"="mds" , "NMDS"="nmds", 
                                                 "PCoA"="pcoa", "CAP"="cap"), selected = "pcoa"),
                         uiOutput("beta.div.cap.param")
                         ),
                  column(width = 3,
                         selectInput("permut.permanova", "Number of permutation", width = "80%", 
                                     choices = seq(100, 1000, by=100)),
                         div(style="display: inline-block;vertical-align:top; width: 50px;", disabled(actionButton("beta.go", "Go!"))),
                         div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("beta.show.dist.tab"))
                         )
                  )
                ),
              fluidRow(
                column(width = 7, offset = 0, style='padding:0px;',
                       box(width = 12, title = "Beta diversity graph", collapsible = T, collapsed = F, solidHeader = T, status = "primary",
                           plotOutput("beta.div.plot"),
                           uiOutput("beta.div.plot.param")
                           )
                       ),
                column(width = 5, offset = 0, style='padding:0px;',
                       box(width = 12, 
                           title = "Pairwise Comparison", collapsible = T, collapsed = F, solidHeader = T, status = "primary",
                           plotOutput("beta.div.pair.plot"),
                           uiOutput("beta.div.stat1.param")
                       )
                )
              )
              ),
      tabItem(tabName = "rb_tab", 
              
                fluidRow(
                  column(width = 12, offset = 0, style='padding:0px;',
                         box(
                           width = 12, title = "Stacked barchart", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                             fluidRow(
                               column(width = 4,
                                      selectInput("st.barchart.factor", "Factor", choices = "choose"),
                                      checkboxInput("st.barchart.change.order", label = "Change X-axis order", value = F),
                                      uiOutput("st.barchart.order")
                                      ),
                               column(width = 4,
                                      selectInput("st.barchart.level.taxa", "Taxa Level", 
                                                  choices = c("Kingdom"="Kingdom", "Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus"),
                                                  selected =  "Class", multiple=TRUE, selectize=TRUE)
                                      ), 
                               column(width = 4,
                                      sliderInput("st.barchart.threshold", "Threashold (%)",min = 0, max = 10, value = 1, step = 0.5),
                                      div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("st.barchart.go", "Go!")),
                                      div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("st.barchart.show.tab"))
                                      )
                               ),
                             plotlyOutput("st.barchart"),
                             uiOutput("param.st.barchart")
                             )
                         )),
              fluidRow(
                column(width = 12, offset = 0, style='padding:0px;',
                       box(width = 12,
                           title = "Heatmap", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                           fluidRow(
                             column(width = 4,
                                    selectInput("heatmap.factor", "Factor", choices = "choose"),
                                    checkboxInput("heatmap.change.order", label = "Change X-axis order", value = F),
                                    uiOutput("heatmap.order")
                             ),
                             column(width = 4,
                                    selectInput("heatmap.level.taxa", "Taxa Level", 
                                                choices = c("Kingdom"="Kingdom", "Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus"),
                                                selected =  "Genus", multiple=TRUE, selectize=TRUE),
                                    checkboxInput("heatmap.rot", label = "Rotate heatmap", value = F)
                             ),
                             column(width = 4,
                                    sliderInput("heatmap.threshold", "Threashold",min = 0, max = 10, value = 1, step = 0.5),
                                    div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("heatmap.go", "Go!")),
                                    div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("heatmap.show.tab"))
                             )
                           ),
                           plotOutput("heatmap", brush = "ht_brush", click = "ht_click"),
                           uiOutput("param.heatmap")
                       )
                )
              ),
              fluidRow(
                  column(width = 12,  offset = 0, style='padding:0px;',
                         box(width = 12, 
                             title = "Barchart", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                             fluidRow(
                               column(width = 4,
                                      selectInput("barchart.factor", "Factor", choices = "choose"),
                                      checkboxInput("barchart.change.order", label = "Change X-axis order", value = F),
                                      uiOutput("barchart.order")
                               ),
                               column(width = 4,
                                      selectInput("barchart.level.taxa", "Taxa Level", 
                                                  choices = c("Kingdom"="Kingdom", "Phylum"="Phylum", "Class"="Class", "Order"= "Order", "Family"="Family", "Genus"="Genus"),
                                                  selected =  "Genus", multiple=TRUE, selectize=TRUE)
                               ),
                               column(width = 4,
                                      sliderInput("barchart.num.OTUs", "Number of OTUs",min = 0, max = 15, value = 1, step = 1),
                                      div(style="display: inline-block;vertical-align:top; width: 50px;",actionButton("barchart.go", "Go!")),
                                      div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("barchart.show.tab"))
                                      
                               )
                             ),
                             plotlyOutput("barchart"),
                             uiOutput("param.barchart")
                             )
                         )
                  )
              ),
      tabItem(tabName = "alpha_tab",
              box(
                width = 12, title = "Alpha diversity parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                fluidRow(
                  column(width = 4, 
                         selectInput("alpha.div.factor", "Factor", choices = c("choose..")),
                         checkboxInput("alpha.div.change.order", label = "Change X-axis order", value = F),
                         uiOutput("alpha.div.order")
                         ),
                  column(width = 4, 
                         selectInput("alpha.div.indices", "Alpha diversity indice(s)", 
                                     choices = c("Goods coverage" = "goods", "Observed OTUs" = "obs.otus",
                                                 "ACE" = "ace", "Chao" = "chao", 
                                                 "Simpson" = "simpson", "Shannon" = "shannon", 
                                                 "Evenness" = "evenness", 
                                                 "PD"= "pd"), 
                                     selected = "obs.otus", multiple=TRUE, selectize=TRUE)
                         ),
                  column(width = 4, 
                         numericInput("alpha.pval.cutoff", "P-value cutoff",value =  0.05, min = 0.02, max = 1, step = 0.01),
                         div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("alpha.div.go", "Go!")),
                         div(style="display: inline-block;vertical-align:top; width: 50px;", uiOutput("alpha.div.show.tab"))
                         )
                  )
                ),
              box(width = 12, boxToolSize = "sm",
                  title = "Boxplot", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                  fluidRow(
                    column(12, align="center", 
                           plotlyOutput("alpha.div.boxplot", width = "80%")
                           )
                    ),
                  uiOutput("param.alpha.div.boxplot")
              
              )
      ),
      tabItem(tabName = "core_tab",
              box(
                width = 12, title = "Core microbiome parameters", collapsible = F, collapsed = F, solidHeader = T, status = "primary",
                fluidRow(
                  column(width = 4, 
                         selectInput("core.factor", "Factor", choices = c("choose.."))
                         
                  ),
                  column(width = 4, 
                         numericInput("core.prevalence", "Prevalence threshold (%)",value =  75, min = 10, max = 100, step = 5)
                  ),
                  column(width = 4, 
                         numericInput("core.rb.detection", "Detection threshold (%)",value =  0.01, min = 0.01, max = 100, step = 0.01),
                         div(style="display: inline-block;vertical-align:top; width: 50px;", actionButton("core.go", "Go!"))
                  )
                )
              ),
              fluidRow(
                column(width = 7, offset = 0, style='padding:0px;',
                       box(width = 12, boxToolSize = "sm",
                           title = "Presence/Absence table", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                  fluidRow(
                    column(12, align="center", 
                           div(style = 'overflow-x: scroll',DT::dataTableOutput("core.table",width = "100%"))
                    )
                  ),
                  uiOutput("param.core.table")
                  )),
                column(width = 5, offset = 0, style='padding:0px;',
                       box(width = 12, boxToolSize = "sm",
                           title = "Venn Diagram", collapsible = T, collapsed = T, solidHeader = T, status = "primary",
                           fluidRow(
                             column(12, 
                                    plotOutput("core.venn.plot"),
                                    uiOutput("core.venn.param")
                             )
                             )
                           )
                       ))
              
      )
    )
)

)

