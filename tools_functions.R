
#########################################
# create taxonomy matrix
#########################################
split.taxon <- function(taxamat){
  list_taxon <- as.vector(taxamat$Taxon)
  split_taxon <- strsplit(list_taxon, ";") 
  taxa_tab <- data.frame(do.call(rbind.data.frame, split_taxon))
  taxa_tab <- taxa_tab[,1:6]
  colnames(taxa_tab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  rownames(taxa_tab) <- taxamat[,1]
  
  return(taxa_tab)
}

#########################################
# Select factors
#########################################
change.factor.order <- function(metamat, factor){
  dat <- as.factor(metamat[, factor])
  factor_var <- levels(dat)
  return(factor_var)
}


#########################################
# Select factors based on physeq object -alpha div only-
#########################################
change.factor.order2 <- function(physeq, factor.tab){
  metamat <- sample_data(physeq)
  dat <- as.factor(metamat[[factor.tab]])
  factor_var <- levels(dat)
  return(factor_var)
}

#######################################
# Functions for sorting factor levels #
# (janhove.github.io)                 #
#######################################

# Sort factor levels by the factor level mean of another covariate
sortLvlsByVar.fnc <- function(oldFactor, sortingVariable, ascending = TRUE) {
  
  require("dplyr")
  require("magrittr")
  
  # Combine into data frame
  df <- data.frame(oldFactor, sortingVariable)
  
  ###
  ### If you want to sort the levels by, say, the median, sd etc. instead of the mean,
  ### just change 'mean(sortingVariable)' below to, say, 'median(sortingVariable)'.
  ###
  
  # Compute average of sortingVariable and arrange (ascending)
  if (ascending == TRUE) {
    df_av <- df %>% group_by(oldFactor) %>% summarise(meanSortingVariable = mean(sortingVariable)) %>% 
      arrange(meanSortingVariable)
  }
  
  # Compute average of sortingVariable and arrange (descending)
  if (ascending == FALSE) {
    df_av <- df %>% group_by(oldFactor) %>% summarise(meanSortingVariable = mean(sortingVariable)) %>% 
      arrange(desc(meanSortingVariable))
  }
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = df_av$oldFactor)
  return(newFactor)
}

# Sort factor levels by their frequency of occurrence
sortLvlsByN.fnc <- function(oldFactor, ascending = TRUE) {
  
  require("magrittr")
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = table(oldFactor)  %>% sort(., decreasing = !ascending)  %>% names())
  return(newFactor)
}

# Sort factor levels arbitrarily
sortLvls.fnc <- function(oldFactor, levelOrder) {
  if(!is.factor(oldFactor)) stop("The variable you want to reorder isn't a factor.")
  
  if(!is.numeric(levelOrder)) stop("'order' should be a numeric vector.")
  
  if(max(levelOrder) > length(levels(oldFactor))) stop("The largest number in 'order' can't be larger than the number of levels in the factor.")
  
  if(length(levelOrder) > length(levels(oldFactor))) stop("You can't have more elements in 'order' than there are levels in the factor.")
  
  if(length(levelOrder) == length(levels(oldFactor))) {
    reorderedFactor <- factor(oldFactor, levels = levels(oldFactor)[levelOrder])
  }
  
  if(length(levelOrder) < length(levels(oldFactor))) {
    levelOrderAll <- c(levelOrder, (1:length(levels(oldFactor)))[-levelOrder])
    reorderedFactor <- factor(oldFactor, levels = levels(oldFactor)[levelOrderAll])
  }
  
  return(reorderedFactor)
}

#########################################
# order factors in meta data
#########################################

order.levels.df <- function(df, factor, orderd_factor1, changefac){
  if (changefac == T){
    Ids <- as.factor(df[,c(factor)])
    level_position <- match(orderd_factor1, levels(Ids))
    Ids <- sortLvls.fnc(Ids, level_position)
    df2<- cbind(df[, -which(names(df) %in% factor)], Ids)
    colnames(df2) <- colnames(df) 
    return(as.data.frame(df2))}
  else return(df)
  
}

###############################################
# Combine Relative abundance and Standard error
###############################################

merge.tab.se.fun <- function(tab, se.tab, add.se){
  
  if(add.se == T){
    var <- colnames(tab) 
    #se.tab <- se.tab[,var]
    merged.tab <-list()
    for(i in 1:length(var)){
      merged.tab[[i]] <- paste(tab[,var[i]], se.tab[,var[i]], sep = "\u00B1")
      }
    merged.tab <- as.data.frame(merged.tab, row.names = rownames(tab))
    colnames(merged.tab) <- colnames(tab)
    tab.out <- merged.tab
  }
  else
    tab.out <- tab
  
  return(tab.out)  
}
