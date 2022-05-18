# Qplots Version 3.0
# Last modified on 2/05/2019
# Script: Guide (server functions)
# Author: Naima BELMOKHTAR (Email: naima1503@gmail.com)
#*******************************************************

#########################################
# Guide

output$home <- renderUI({includeHTML("www/home.html")})
output$mgem <- renderUI({includeHTML("www/mgem.html")})

output$qplot_workflow <- renderUI({
  bsCollapse(id = "collapseqplot", open = "",
             bsCollapsePanel("Overview", htmlOutput("overview_pipeline"), style = "primary"),
             bsCollapsePanel("Import Data", htmlOutput("import_data"), style = "primary"),
             bsCollapsePanel("Beta Diversity", htmlOutput("betaguide"),style = "success"),
             bsCollapsePanel("Alpha Diversity", htmlOutput("alphaguide"),style = "warning"),
             bsCollapsePanel("Relative Abundance", htmlOutput("rb_guide"), style = "danger"),
             bsCollapsePanel("Core microbiome", htmlOutput("rb_guide"), style = "info")
             
             #,
             #bsCollapsePanel("Serial group comparison", htmlOutput("serialguide"),style = "info"),
             #bsCollapsePanel("Network Input", htmlOutput("netinputguide"),style = "info"),
             #bsCollapsePanel("Co-occurrence Network", htmlOutput("overview_net"),
              #               bsCollapsePanel("Cytoscape", htmlOutput("constract_net") , style = "primary"),
               #              bsCollapsePanel("Gephi", htmlOutput("vis_gephi") , style = "primary"), 
                #             style = "primary")
             )
})

output$overview_pipeline <- renderUI(includeHTML("www/overview_workflow.html"))
output$import_data <- renderUI(includeHTML("www/data_input.html"))
output$rb_guide <- renderUI(includeHTML("www/r_abundance.html"))
output$alphaguide <- renderUI(includeHTML("www/alpha_div.html"))
output$betaguide <- renderUI(includeHTML("www/beta_div.html"))
output$serialguide <- renderUI(includeHTML("www/serial_comp.html"))
output$netinputguide <- renderUI(includeHTML("www/input_network.html"))
output$constract_net <- renderUI(includeHTML("www/constract_net.html"))
output$vis_gephi <- renderUI(includeHTML("www/gephi_net.html"))
output$overview_net <- renderUI(includeHTML("www/requirement_net.html"))

#Fastq Analysis Guide
output$script_guide <- renderUI({
  bsCollapse(id = "collapsefq", open = "",
             bsCollapsePanel("Overview", htmlOutput("overview_fqanalysis"), style = "primary"),
             bsCollapsePanel("Requirments", htmlOutput("fq_req"), style = "info"),
             bsCollapsePanel("Useful commands", htmlOutput("fq_u_command"), style = "success"),
             bsCollapsePanel("File Transfer", htmlOutput("fq_trans"), style = "danger")
             
  )
})

output$overview_fqanalysis <- renderUI(includeHTML("www/overview_fq.html"))
output$fq_req <- renderUI(includeHTML("www/fq_requirment.html"))
output$fq_u_command <- renderUI(includeHTML("www/u_command.html"))
output$fq_trans <- renderUI(includeHTML("www/fq_transfer.html"))

#Alpha diversity Qiime Guide

output$alpha_guide <- renderUI({
  bsCollapse(id = "collapse2", open = "", multiple = F,
             bsCollapsePanel("Overview", htmlOutput("overview_alpha") , style = "success"),
             bsCollapsePanel("Input Files", htmlOutput("inputfiles_alpha") , style = "info"),
             bsCollapsePanel("Hands-On", htmlOutput("handson_alpha"),
                             bsCollapsePanel("Required Files", htmlOutput("requiredfile_alpha") , style = "primary"),
                             bsCollapsePanel("Upload Files / arguments", htmlOutput("uploadfile_alpha") , style = "primary"),
                             bsCollapsePanel("Results", htmlOutput("results_alpha") , style = "primary"), 
                             style = "primary")
  )
})

output$overview_alpha <- renderUI(includeHTML("www/overview_alpha.html"))
output$inputfiles_alpha <- renderUI(includeHTML("www/inputfiles_alpha.html"))
output$handson_alpha <- renderUI(includeHTML("www/handson_alpha.html"))
output$requiredfile_alpha <- renderUI(includeHTML("www/requiredfile_alpha.html"))
output$uploadfile_alpha <- renderUI(includeHTML("www/uploadfile_alpha.html"))
output$results_alpha <- renderUI(includeHTML("www/results_alpha.html"))

