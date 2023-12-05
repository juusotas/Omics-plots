#### Functions to plot gene set analysis results ####
### TO DO ####
### Paralel plot_gsea_results loop for faster execution ###
### gsea_volcano needs an option to plot only one of the plots, i.e. if only positive or negative NES is present

### Function to mine Log2 Fold Change information which is attached to plots using other functions ###
get_gene_information = function(result_df,gene_df,symbol_col, FC_col){
  # result_df = gene set analysis result
  # gene_df = data frame with symbol and log2 fold change
  # symbol_col = col index of symbol col in gene_df
  # FC_col = col index of log2fold change col in gene_df
  out_df = data.frame()
  gene_df = cbind(gene_df[symbol_col],gene_df[FC_col])
  if("core_enrichment" %in% colnames(result_df)){
    genes = result_df$core_enrichment
  }else{
    genes = result_df$geneID
  }
  genes = gsub("/", ";",genes) 
  genes = strsplit(genes, ";")
  for(i in 1:(nrow(result_df))){
    gene_set_FC=gene_df[gene_df[[1]] %in% genes[[i]],]
    gene_set_rep = unlist(result_df[i,])
    gene_set_inf = data.frame(matrix(rep(gene_set_rep, each = nrow(gene_set_FC)), nrow(gene_set_FC)))
    temp_df = cbind(gene_set_inf,gene_set_FC)
    out_df = rbind(out_df,temp_df)
  }
  
  name_cols = c(colnames(result_df),"symbol","log2fc")
  colnames(out_df) = name_cols
  return(out_df)
}

### Plotting 
ridge_plot = function(geneset_df, n_pathways){
  # geneset_df = resulting data frame from get_gene_info function,
  # n_pathways = how many pathways to plot
  asd_t = unique(geneset_df$Description)
  asd_t = asd_t[1:n_pathways]
  asd_t1 = geneset_df[geneset_df$Description %in% asd_t,]
  
  
  gg1 = ggplot(asd_t1, aes(x = log2fc, y = Description,fill=padjust))+ 
    geom_density_ridges(
      jittered_points = TRUE, scale = .95, rel_min_height = .01,
      point_shape = "*", point_size = 3, size = 0.25,
      position = position_points_jitter(height = 0)) +
    labs(title = 'Log2 Fold Change distribution of genesets', x = "Log2 Fold Change", y = "Gene sets") + 
    scale_fill_gradient(high = "blue", low = "red", name = "Adjusted p-value (BH)")
  return(gg1)
}

violin_plotly = function(geneset_df,n_pathways){
  # geneset_df = resulting data frame from get_gene_info function,
  # n_pathways = how many pathways to plot
  asd_t = unique(geneset_df$Description)
  asd_t = asd_t[1:n_pathways]
  df = geneset_df[geneset_df$Description %in% asd_t,]
  
  fig = df %>% 
    plot_ly(x = ~Description, y = ~log2fc, split =~Description, 
            hovertext = ~symbol, type = 'violin', points='all',hoveron= "points", 
            showlegend=FALSE, box = list(visible = T), meanline = list(visible = T)
    )
  return(fig)
}

violin_plot = function(geneset_df, n_pathways){
  # geneset_df = resulting data frame from get_gene_info function,
  # n_pathways = how many pathways to plot
  asd_t = unique(geneset_df$Description)
  asd_t = asd_t[1:n_pathways]
  df = geneset_df[geneset_df$Description %in% asd_t,]
  
  p <- ggplot(df, aes(y=Description, x=log2fc,fill=Description)) + 
    geom_violin(trim=F) + theme(legend.position="none") +
    geom_boxplot(width=0.1, fill="white")+
    labs(x = "Log2 Fold Change",y = "Gene Sets")
  p
  
}
gsea_volcano = function(gsea_result,n_pathways){
  # gsea_result = a gsea_result table from for example fgsea function
  # n_pathways = how many pathways to plot
  # add functionality for only one plot
  x = gsea_result
  x["regulation"] = ifelse(x$NES > 0,"positive","negative")
  x_up = x[x$regulation == "positive",]
  x_up = mutate(x_up,Description = fct_reorder(Description, NES))
  x_down = x[x$regulation == "negative",]
  x_down = mutate(x_down,Description = fct_reorder(Description, NES,dplyr::desc))
  
  if(nrow(x_up) > n_pathways){
    x_up = x_up[1:n_pathways,]
  }
  if(nrow(x_down) > n_pathways){
    x_down= x_down[1:n_pathways,]
  }
  
  nes_up = ggplot(x_up, aes(x=NES, y=Description)) +
    geom_point(aes(colour = p.adjust),size=4, alpha=0.7)+
    scale_y_discrete(position = "right")+
    theme(axis.title.y = element_blank())+
    scale_colour_gradient2(
      mid = "blue",
      high = "red",
      space = "Lab",
      guide = "colourbar",
      aesthetics = "colour"
    )
  nes_down = ggplot(x_down, aes(x=NES, y=Description)) +
    geom_point(aes(colour = p.adjust),size=4, alpha=0.7)+
    scale_y_discrete(position = "left")+
    theme(axis.title.y = element_blank())+
    scale_colour_gradient2(
      mid = "blue",
      high = "red",
      space = "Lab",
      guide = "colourbar",
      aesthetics = "colour"
    )
  patch = ggarrange(nes_down,nes_up, common.legend =T,legend = "bottom")
  return(patch)
}

bar_plot = function(geneDF, n_categories){
  geneset_df = geneDF@result
  if("Count" %in% colnames(geneset_df)){
    set_size = geneset_df$Count
  }else{
    set_size = geneset_df$setSize
  }
  gene_sub_set = geneset_df[1:n_categories,]
  set_size = set_size[1:n_categories]
  gg = ggplot(data=gene_sub_set, aes(x=set_size, y=Description, fill=p.adjust)) +
    geom_bar(stat="identity") + theme(text = element_text(size=10))
  return(gg)
}

#### This fucntion is used to loop through a list of gene set analysis results, which are then as inputs for different plot functions ####
#### Such as thoes above ####
plot_gsea_results = function(gene_data,analysis_list,plot_anno_list,l2fc = "yes", npathways = 30,n_ema = 50,padj=0.05){
  # gene_data = data frame with log2foldchange and symbol columns
  # analysis_list = list of gene set analysis results
  # plot_anno_list = list of annotation strings for each plot, for exmaple the plot title or plot file name
  # l2fc = if "yes" make plots that use log2 fold change information, if "no" make plots without log2 fold change information
  # npathways = How many pathways to plot
  # n_ema = How many pathways to plot in the emaplot function
  # padj = What is the padjusted threshold for some plots
  library(forcats)
  for(z in 1:length(analysis_list)){
    options(ggrepel.max.overlaps = Inf)
    y = analysis_list[[z]]
    print(nrow(y@result))
    if(nrow(y@result) > 0 ){
      if(l2fc == "no"){
        plot_title = paste("Top",n,plot_anno_list$plotheader[z])
        
        bar.plot = bar_plot(y,n_categories =  n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10))
        ggsave(plot_anno_list$barplotfilename[z],plot = bar.plot, scale = 1.5, width = 24,height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        dot_plot = dotplot(y, x= "GeneRatio", color = "p.adjust",showCategory = n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10))
        ggsave(plot_anno_list$dotplotfilename[z],plot = dot_plot, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        upset.plot = upsetplot(y,n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10))
        ggsave(plot_anno_list$upsetplotfilename[z],plot = upset.plot, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
        
      }
      if(l2fc == "yes"){
        cnet_title  = paste("Top",10,plot_anno_list$cnetheader[z])
        gg = cnetplot(y, categorySize="pvalue", color.params = list(foldChange = GSEA),showCategory = 10, node_label = "all",layout = "dh")+ 
          labs(title= cnet_title)+theme(plot.title = element_text(hjust = 0.5, size = 10))
        ggsave(plot_anno_list$cnetfilename[z],plot= gg, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 200, device = "tiff")

        y_result = subset(y@result,p.adjust <= padj)
        gene_info = get_gene_information(y_result,gene_data,1,2)
        colnames(gene_info)[6] = "padjust"
        gene_info = transform(gene_info,padjust = as.numeric(padjust))
        row_gene_info = nrow(y_result)
        
        if( row_gene_info < npathways ){
          n = row_gene_info
        }else{
          n = npathways
        }
        
        plot_title = paste("Top",n,plot_anno_list$plotheader[z])
        
        
        fig = violin_plotly(gene_info,n)
        fig = fig %>% layout(title = plot_title, xaxis = list(title = "Gene set"), yaxis = list(title ="Log2 Fold Change"))
        htmlwidgets::saveWidget(as_widget(fig), plot_anno_list$violinplotfilename[z])
        
        viol_plt = violin_plot(gene_info,n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10 ))
        
        filename = plot_anno_list$violinplotfilename[z]
        filename = str_replace(filename, ".html", ".tiff")
        ggsave(filename,plot = viol_plt, scale = 1, width = 26, height =16 ,units = "cm",dpi = 600, device = "tiff")
        
        
        ridge = ridge_plot(gene_info,n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10 ))
        ggsave(plot_anno_list$ridgeplotfilename[z],plot = ridge, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        bar.plot = bar_plot(y,n_categories =  n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10))
        ggsave(plot_anno_list$barplotfilename[z],plot = bar.plot, scale = 1.5, width = 24,height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        dot_plot = dotplot(y, x= "GeneRatio", color = "p.adjust",showCategory = n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10))
        ggsave(plot_anno_list$dotplotfilename[z],plot = dot_plot, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        upset.plot = upsetplot(y,n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10))
        ggsave(plot_anno_list$upsetplotfilename[z],plot = upset.plot, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        
        heat = heatplot(y,foldChange = GSEA,showCategory = n)+
          labs(title= plot_title)+theme(plot.title = element_text(hjust = 0.5, size = 14),
                                        axis.text.x =element_text(size = 10, angle = 90 ))
        ggsave(plot_anno_list$heatmapfilename[z],plot= heat, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
        
        
        
        tryCatch({
          
          y_pw_termsim = pairwise_termsim(y)
          ema = emapplot(y_pw_termsim,showCategory = n_ema , layout.params = list(layout = "kk", coords = NULL), color = "p.adjust", 
                         cex.params = list(category_node = 0.5, category_label = 0.5, line = 0.3),
                         hilight.params = list(catergory = NULL,alpha_hilight = 1, alpha_no_hilight = 0.3),
                         cluster.params = list(cluster = T, method = stats::kmeans, n = NULL, legend =
                                                 T, label_style = "shadowtext", label_words_n = 4, label_format = 30))+
            labs(title= plot_title)+theme(plot.title = element_text(hjust = 1, size = 15))
          ggsave(plot_anno_list$emplotfilename[z],plot= ema, scale = 1.5, width = 24, height =24 ,units = "cm",dpi = 600, device = "tiff")
          
          treefilename = str_replace(plot_anno_list$emplotfilename[z],"emaplot","treeplot")
          
          tree = treeplot(y_pw_termsim,showCategory = n,color = "p.adjust", cluster.params = list(method = "ward.D2", n = 5))+
            labs(title= plot_title)+theme(plot.title = element_text(hjust = 1, size = 15))
          ggsave(treefilename,plot= tree, scale = 1.5, width = 24, height =12 ,units = "cm",dpi = 600, device = "tiff")
          
        },
        warning=function(war)
        {
          print(paste("MY WARNING: ", war))
        },
        error=function(err)
        {
          print(paste("MY ERROR: ", err))
        },
        finally=function(f)
        {
          print(paste("e: ", e))
        })

        if(class(y)[1] == "gseaResult"){
          tryCatch({
            gene_sub_set = y@result
            NES_plot_title = paste("Top",n,"NES of",plot_anno_list$plotheader[z])
            NES_file_name = paste("NES_barplot",plot_anno_list$plotheader[z],".tiff")
            vol_file_name = paste("gseaVolcano",plot_anno_list$plotheader[z],".tiff")
            
            gene_sub_set$Description = substr(as.matrix(gene_sub_set$Description), 1,100)
            colnames(gene_sub_set)[2] = "Description"
            gseavol = gsea_volcano(gene_sub_set,60)
            gseavol = annotate_figure(gseavol, top = text_grob(plot_anno_list$plotheader[z],color = "black", size = 14))
            ggsave(vol_file_name,plot= gseavol, scale = 1.5, width = 24, height =12,units = "cm",dpi = 200, device = "tiff")
            
            gene_sub_set = gene_sub_set %>% arrange(desc(abs(NES)))
            gene_sub_set = slice_head(gene_sub_set, n= n)
            # Reorder following the value of another column:
            gg_nes = gene_sub_set %>%
              mutate(Description = fct_reorder(Description, NES)) %>% ggplot(aes(x=NES, y=Description, fill=p.adjust)) +
              geom_bar(stat="identity") + theme(text = element_text(size=15)) +
              labs(title=NES_plot_title )+theme(plot.title = element_text(hjust = 0.5, size = 15),
                                                axis.text.x =element_text(size = 10, angle = 90 ))
            ggsave(NES_file_name,plot= gg_nes, scale = 1.5, width = 24, height =12,units = "cm",dpi = 200, device = "tiff")
            
            fig_nes = violin_plotly(gene_info,n)
            fig_nes = fig_nes %>% layout(title = paste(plot_title,"NES",sep = " "), xaxis = list(title = "Gene set"), yaxis = list(title ="Log2 Fold Change"))
            htmlwidgets::saveWidget(as_widget(fig_nes), paste("NES",plot_anno_list$violinplotfilename[z],""))
          },
          warning=function(war)
          {
            print(paste("MY WARNING: ", war))
          },
          error=function(err)
          {
            print(paste("MY ERROR: ", err))
          },
          finally=function(f)
          {
            print(paste("e: ", e))
          })
          
        }
        
      }
    }  
    else{
      next
    }
  }
}
