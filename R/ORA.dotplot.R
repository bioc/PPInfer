# ORA dotplot
ORA.dotplot <- function(object, category, size, count, pvalue, top=10,
                        sort = NULL, decreasing = FALSE,
                        p.adjust.methods = NULL, numChar = NULL,
                        title = NULL, transparency = 0.5, plot = TRUE)
{
  object <- data.frame(object)
  
  if(!is.null(sort))
  {
    index <- order(object[,sort], decreasing = decreasing)
    object <- object[index,]
  }
  
  EnrichTab <- object[1:min(top,nrow(object)),]
  EnrichTab <- EnrichTab[,c(category, size, count, pvalue)]
  colnames(EnrichTab) <- c('Category', 'Size', 'Count', 'Pvalue')
  
  if(is.null(numChar))
  {
    numChar <- max(nchar(as.character(EnrichTab$Category)))
  }
  else
  {
    if(length(unique(substr(EnrichTab$Category, 1, numChar))) < nrow(EnrichTab))
    {
      numChar <- max(nchar(as.character(EnrichTab$Category)))
      message('Note : numChar is too small.', '\n')
    }
  }
  EnrichTab$Category <- paste(substr(EnrichTab$Category, 1, numChar),
                           ifelse(nchar(as.character(EnrichTab$Category)) > numChar,
                                  "...", ""), sep = "")
  GeneRatio <- EnrichTab$Count/EnrichTab$Size
  
  index <- which(colnames(object) == pvalue)
  
  if(!is.null(p.adjust.methods))
  {
    adj.Pvalue <- p.adjust(object[,index], p.adjust.methods)
    adj.Pvalue <- adj.Pvalue[1:nrow(EnrichTab)]
    EnrichTab <- data.frame(EnrichTab, GeneRatio, adj.Pvalue)
  }
  else
  {
    EnrichTab <- data.frame(EnrichTab, GeneRatio)
  }
  
  EnrichTab$Category <- factor(EnrichTab$Category,
                               levels = EnrichTab$Category[nrow(EnrichTab):1])
  
  p <- ggplot(EnrichTab, aes_string(x='GeneRatio', y='Category', size='Count', 
                                 color=ifelse(!is.null(p.adjust.methods) &
                                                is.numeric(object[,index]),
                                              'adj.Pvalue', 'Pvalue'))) + 
    geom_point(alpha = transparency) + ylab('')
  
  if(plot == TRUE)
  {
    print(EnrichTab); return(p)
  }
  else
  {
    print(p); return(EnrichTab)
  }
}