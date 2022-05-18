# Qplots Version 1.0
# Script: Beta Diversity (Analysis functions)
# Author: Naima BELMOKHTAR (Email: naima1503@gmail.com)
#*******************************************************

#########################################
# import requested functions
source("tools_functions.R")

#Calculate beta diversity (output: Distance matrix )
beta.dist <- function(physeq, dist.method){
  #input data
  otu_file <- otu_table(physeq)
  otu_file <- as.data.frame(otu_file@.Data)
  otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
  
  tree_file <- phy_tree(physeq)
  
  #prepare data
  otu_file <- otu_file[,order(names(otu_file))]
  otu_file <- data.frame(t(otu_file))
  
  #Calculate distance matrix based on Unifract - Bray-curtis
  # Bray-curtis
  if(dist.method=="bray")
    Dist_matrix <- as.matrix(vegdist(otu_file, method="bray", diag = T, upper = T))
  # Unifract including: Weighted, unweighted, Gunifract
  else{
    unifracs <- GUniFrac(otu_file, tree_file, alpha = c(0.0,0.5,1.0))$unifracs
    Dist_matrix <- unifracs[, , dist.method]
  }
  
  return(Dist_matrix)
}

#generate all group (output: factor of group per sample)
all.groups <- function(physeq, group_name){
  #input data
  meta_file <- sample_data(physeq)
  mappingVar <- names(meta_file)
  meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
  
  #prepare data
  meta_file <- data.frame(meta_file[order(row.names(meta_file)),])
  names(meta_file) <- mappingVar
  meta_file_pos <- which(colnames(meta_file) == group_name)
  
  all_groups <- as.factor(meta_file[,meta_file_pos])
  return(all_groups)
}

#Ordination plots: MDS, NMDS, PCoA, CAP
beta.ordi <- function(physeq, dist.method="d_0.5", div.factor, ordination.method="cap", m=11, permut=100){
  
  dis.matrix <- beta.dist(physeq, dist.method)
  
  all_groups <- all.groups(physeq, div.factor)
  
  dist_comp <- dis.matrix[, !is.na(all_groups)]
  
  
  meta_file <- sample_data(physeq)
  meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
  meta_file <- data.frame(meta_file[order(row.names(meta_file)),], stringsAsFactors = T)
  meta_file[, div.factor] <- as.factor(meta_file[, div.factor])
  
  if(length(unique(meta_file[, div.factor]))>=2){
   
    adonis.res <-adonis2(as.formula(paste("as.dist(dist_comp) ~", div.factor, sep = "")), data = meta_file)
     if(ordination.method == "pcoa"){
       Ordi_pcoa <- pcoa(dist_comp, correction="none", rn=NULL)
       colnames(Ordi_pcoa[["vectors"]]) <- paste0("Dim",1:dim(Ordi_pcoa[["vectors"]])[2])
       
       # Dimension 
       Axis1.percent <- Ordi_pcoa$values$Relative_eig[[1]] * 100
       Axis2.percent <- Ordi_pcoa$values$Relative_eig[[2]] * 100
       
       plot1 <- ordiplot(Ordi_pcoa$vectors, choices = c(1,2), display = "sites")
       xlabel = paste("PCo1 (", round(Axis1.percent, 2), " %", "of total variation)", sep = "") 
       ylabel = paste("PCo2 (", round(Axis2.percent, 2), " %", "of total variation)", sep = "")
       } else if(ordination.method == "cap"){
         #meta_file[,"Generation"] <- factor(meta_file[,"Generation"])
         #Ordi_cap <- CAPdiscrim2(formula = as.formula(paste("as.dist(dist_comp) ~", div.factor)), data = meta_file,  m = 11, axes = 2, permutations =permut)
         Ordi_cap <- CAPdiscrim2(x = as.dist(dist_comp), group = div.factor , data = meta_file,  m = m, axes = 2, permutations =as.numeric(permut))
    
         plot1 <- ordiplot(Ordi_cap, display = "sites")
         xlabel <- "CAP1"
         ylabel <- "CAP2"
         } else if(ordination.method == "nmds") {
           # Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
           Ordi_nmds <- metaMDS(dist_comp, k = 2, trymax = 100)
           
           plot1 <- ordiplot(Ordi_nmds, display = "sites")
           xlabel <- "nMDS1"
           ylabel <- "nMDS2"
           } else if(ordination.method == "mds"){
             # Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
             Ordi_mds <-cmdscale(dist_comp, k = 2)
             plot1 <- ordiplot(Ordi_mds, display = "sites")
             xlabel <- "nMDS1"
             ylabel <- "nMDS2"
             }
    sites.long1 <- data.frame("Factor1"=meta_file[[div.factor]], "axis1"=plot1$sites[, 1], "axis2"=plot1$sites[, 2], "labels"=rownames(plot1$sites))
    colnames(sites.long1)[1]<- paste(div.factor)
    beta.ordi.out <- list("ordi.tab" = sites.long1,
                          "xlabel" = xlabel,
                          "ylabel" = ylabel,
                          "meta_file.factor" = meta_file[, div.factor],
                          "adonis" = adonis.res$`Pr(>F)`[1] #adonis.res[[1]][6][[1]][1]
                          )
    return(beta.ordi.out)
  } else return(NULL)
}

beta.pair.basic <- function(physeq, dist.method="d_0.5", div.factor, permut=100){
  
  dis.matrix <- beta.dist(physeq, dist.method)
  
  all_groups <- all.groups(physeq, div.factor)
  
  dist_comp <<- dis.matrix[, !is.na(all_groups)]
  
  
  meta_file <- sample_data(physeq)
  meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),],row.names=row.names(meta_file))
  meta_file <- data.frame(meta_file[order(row.names(meta_file)),], stringsAsFactors = T)
  meta_file[, div.factor] <- as.factor(meta_file[, div.factor])
  
  if(length(unique(meta_file[, div.factor]))>2){ 
    
    pairwise.stat.test <- pairwise.adonis2(as.formula(paste("as.dist(dist_comp) ~", div.factor, sep = "")), data = meta_file)
    
    
    df <- as.data.frame(do.call(cbind.data.frame, pairwise.stat.test))
    df2 <- t(df[1,2:dim(df)[2]])
    df2 <- data.frame("Groups"=rownames(df2), df2, row.names = NULL)
    df2 <- df2 %>% 
      separate(Groups, c("Groups", "values"), sep = "\\.")
    
    df3 <- dcast(df2, Groups ~ values, value.var = paste(div.factor))
    
    
    df.select <- df3 %>% 
      separate(Groups, c("Groups1", "Groups2"), sep = "_vs_") %>%
      select(Groups1, Groups2, "Pval"="Pr(>F)")
    
    #Plot tringle heatmap of p-values
    df.select.dcast <- dcast(df.select, Groups1 ~ Groups2, value.var = "Pval")
    df.select.dcast2 <- data.frame(df.select.dcast[,2:dim(df.select.dcast)[2]], row.names = df.select.dcast$Groups1, check.names = F)
    df.select.dcast2.ord <- df.select.dcast2[order(rowSums(is.na(df.select.dcast2))), order(colSums(is.na(df.select.dcast2))) ]
    
    df.select$Groups1 <- factor(df.select$Groups1, levels = rownames(df.select.dcast2.ord))
    df.select$Groups2 <- factor(df.select$Groups2, levels = colnames(df.select.dcast2.ord))
    df.select$p_range <- cut(df.select$Pval,               # Add group column
                             breaks = c(0, 0.001, 0.01, 0.05, 1),
                             labels = c("p<0.001", "0.001<p<=0.01", "0.01<p<=0.05", "p>0.05"))
    
    #col_breaks <- c("(0,0.001]" = "#046C9A", "(0.001,0.01]" = "#3B9AB2", "(0.01,0.05]" = "#78B7C5", "(0.05,1]" = "#EBCC2A")
    col_breaks <- c("p<0.001" = "#046C9A", "0.001<p<=0.01" = "#3B9AB2", "0.01<p<=0.05" = "#78B7C5", "p>0.05" = "#EBCC2A")
    
    beta.paiwise.stat1.out <- list("df.select" = df.select,
                          "col_breaks" = col_breaks,
                          "dist.p.values" = df.select.dcast2.ord
                          )
    return(beta.paiwise.stat1.out)
  } else
    return(NULL)

  
}
  
