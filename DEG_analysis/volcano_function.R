#args - lab_type, sig #top options

#lab_type = NULL (as in they dont put anything in)
#lab_type = genes
#lab_type = sig_genes
#lab_type = top_genes 
#top = 20, if(top(is.missing()) #only a problem if it's evaluated so if it's missing but we don't use it in the other settings no trouble, default can be modified.


#Example scenarios:

plot_vol(deg, p = 0.01, title = "D28-Postvac") #No lab and no genes
plot_vol(deg, p = 0.01, lab = 'none', title = "D28-Postvac") #Specified no labs
plot_vol(deg, lab = 'sig', genes = c("FCRL", "APOL"), title = "D28-Postvac")
plot_vol(deg, genes = c("FCRL", "APOL"), title = "D28-Postvac") #No lab setting declared, but genes is.character -> use sig as default lab setting

plot_vol(deg, lab_type = 'top', genes = c("FCRL", "APOL"), title = "D28-Postvac") #label all genes specified

plot_vol(deg, p = 0.01, lab = 'top', title = "D28-Postvac") #label top 20 DEGS
plot_vol(deg, p = 0.01, lab = 'top', top = 50, title = "D28-Postvac") #top 50


#Output testing ####
V1 <- toptables$V1_vs_V0

plot_vol(V1, lab = 'top', title = "24h Postvac")

plot_vol <- function(deg, p = 0.05, FC = 0,
                     lab_type= "top", genes, top = 20,
                     gene_col = "gene_name", #colname with gene_IDs
                     title="", alpha = 0.98,
                     colours = c("#a50000","#800000","#ef5a3a","orange","yellow")) {

  
#/ Checks
  #missing(deg) || stop("Error: deg arg is missing. Please provide a toptable data frame")
  #!missing(genes) & !is.character(genes) || stop("Error: label is not a character vector")

#Also stop if lab_type is present but no genes for labs
if (!missing(lab_type) & lab_type != c('top') & missing(genes)) {
  stop("Error: label requires character vector with selected genes")}

  #/ 1) Define sig values
  log_p <- -log10(p)
  #Adj p values #need to convert adjusted p-val to unadjusted for plot yscale
  values <- seq(0.050,0.051, by=0.00001)
  deg$adj.P.Val <- round(deg$adj.P.Val, 5) #This is such a janky way I apologise
  adj.p <- deg[(deg$adj.P.Val %in% (values)),c("P.Value")] 
  adj.p <-  min(adj.p)
  log_adj <- -log10(adj.p)
  
  ##/Set up reg table 2)
 
  deg <- deg %>%
    mutate(reg =
             case_when(
               deg$logFC >= FC & deg$adj.P.Val <= p ~ "Sig Adj. P <0.05",
               deg$logFC <= FC & deg$adj.P.Val <= p ~ "Sig Adj. P <0.05",
               deg$logFC >= FC & deg$P.Value <= p ~ "Sig P <0.05",
               deg$logFC <= FC & deg$P.Value <= p ~ "Sig P <0.05",
               abs(deg$logFC) <= FC & deg$adj.P.Val >= p ~ "No Change",
               abs(deg$logFC) <= FC & deg$adj.P.Val <= p ~ "No Change",
               abs(deg$logFC) > FC & deg$adj.P.Val >p ~ "No Change")) %>%
    mutate(reg =
             factor(reg, levels =
                      c("Sig Adj. P <0.05", "Sig P <0.05", "No Change")))
  
  #Define labels 3)
  
  if (is.null(lab_type)){#No entry for lab list, currently there's no default arg
    gene_label <- c("")
    lab_data <- NULL #Empty dataframe, need to work out how not to error
    warning("No genes highlighted")
  } else if(lab_type == "sig"){
      gene_label <- genes
      lab_data <- deg[(deg$reg == "Sig P <0.05" |deg$reg == "Sig Adj. P <0.05" )
                    & (deg$gene_col %in% genes),]
  } else if (lab_type == "ns"){
        gene_label <- genes
        lab_data <- deg[(deg$gene_col %in% genes),]}
    else if (lab_type == "top"){
          lab_data <- slice_min(deg,adj.P.Val,n=top)
          gene_label <- lab_data$gene_col
  }


#Plot Volcano 4)
vol <-
  deg %>% ggplot(aes(x=logFC,y=-log10(P.Value),label=gene_name))+
  geom_point(aes(color = P.Value, alpha=alpha))+
  labs(title = title) +
  theme_minimal() +theme(legend.position = "none") +
  geom_hline(yintercept = log_p, linetype = 2.5, alpha =0.7) +
  geom_hline(yintercept =log_adj, linetype = 2.5, alpha =0.7)+
  geom_label_repel(data=lab_data,
                   size=3.5,direction="both",nudge_y =1.6,nudge_x =0.1,
                   angle= 70,vjust= 0,segment.size= 0.5,
                   segment.color="#331002",fill="#f7f7f5")+
  scale_color_gradientn(colours = colours,
                        values=c(0,adj.p,p,1))
return(vol)
} 

plot_vol(V1, lab = 'top', title = "24h Postvac")

#Debugging ####
??geom_label_repel



test_data <- slice_min(V1,adj.P.Val,n=20)
test_label <- test_data$gene_name


log_p <- -log10(p)
#Plot volcano
vol <-
  V1 %>% ggplot(aes(x=logFC,y=-log10(P.Value),label=gene_name))+
  geom_point(aes(color = P.Value, alpha=0.98))+
  theme_minimal() +theme(legend.position = "none") +
  geom_hline(yintercept = log_p, linetype = 2.5, alpha =0.7) +
  geom_hline(yintercept =log_adj, linetype = 2.5, alpha =0.7)+
  geom_label_repel(data=test_data, size=3.5,direction="both",nudge_y =1.6,nudge_x =0.1,angle= 70,vjust= 0,segment.size= 0.5,segment.color="#331002",fill="#f7f7f5")+
  scale_color_gradientn(colours = c("#a50000","#800000","#ef5a3a","orange","yellow"),values=c(0,log_adj,log_p,1))
vol
