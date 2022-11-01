

#'Function filters data from missing values and expression lower than certain threshold; Functions works on the output of read functions as exp_data_input and sample_data_input
#' @export


data_filt <- function(){
  #call library
  library(dplyr)


  ######Expression data filtration ###

  exp_data_col <- c("ID_REF", sample_data_input$accession)
  samples_iden_cols <- exp_data_col[-which(exp_data_col =="ID_REF")]

  if (dim(dplyr::filter(exp_data_input,
                        rowSums(exp_data_input[,!colnames(exp_data_input)=="ID_REF"]) > length(samples_iden_cols)))[1]==0){
    exp_data_input <- exp_data_input
  } else{
    #Remove genes with low overall expression
    exp_data_input <- dplyr::filter(exp_data_input, rowSums(exp_data_input[,!colnames(exp_data_input)=="ID_REF"]) > length(samples_iden_cols))
  }


  #Remove missed values
  exp_data_input <- exp_data_input[!exp_data_input$ID_REF=="",]


  #convert ID Ref to NCBI Symbol

  iden_menu_list <- c("ENSGxxx", "ENSPxxx", "Entrez id 'numerical'",
                      "Gene name; e.g.: kinase 1","Probe id; e.g: xxx_at", "Ref seq ID; e.g.: NM_xxx",
                      "Gene symbol; e.g.: DDR1", "IlluminaHT_12_v4 beadchip; e.g: ILMN_xxx ")


  iden_type_list <- c("GENEID", "PROTEINID" , "ENTREZID","GENENAME" ,"PROBEID", "REFSEQ","SYMBOL", "PROBEID")

  print("Please select type of gene identifier submitted")
  print(exp_data_input$ID_REF[1:5])

  iden_type_indx <- menu(iden_menu_list, title= "Select ID_REF type", graphics = TRUE)


  #For loop to extract ensemble ID relevant ID and symbol of NCBI (Entrez database)
  if((iden_type_indx== 5) || (iden_type_list[iden_type_indx]=="REFSEQ") ){
    library(hgu133a.db)
    res_df <- select(hgu133a.db, keys=exp_data_input$ID_REF, columns = c("SYMBOL","ENTREZID"),
                     keytype= iden_type_list[iden_type_indx])
    res_df <- res_df[!duplicated(res_df[,iden_type_list[iden_type_indx]]),]
    exp_data_input$ncbi_id = res_df$ENTREZID
    exp_data_input$ncbi_symbol = res_df$SYMBOL

  } else if (iden_type_indx== 8) {
    library(illuminaHumanv4.db)
    res_df <- AnnotationDbi::select(illuminaHumanv4.db, keys= exp_data_input$ID_REF , keytype = iden_type_list[iden_type_indx],
                                    columns=c("SYMBOL","ENTREZID") )
    res_df <- res_df[!duplicated(res_df[,iden_type_list[iden_type_indx]]),]
    exp_data_input$ncbi_id = res_df$ENTREZID
    exp_data_input$ncbi_symbol = res_df$SYMBOL
  }  else if (iden_type_list[iden_type_indx]=="SYMBOL"){   #end of else if
    library(EnsDb.Hsapiens.v79)
    res_df <- ensembldb::select(EnsDb.Hsapiens.v79, keys= exp_data_input$ID_REF, keytype = "SYMBOL", columns = c("ENTREZID") )
    exp_data_input<- dplyr::left_join(exp_data_input, res_df, by= c("ID_REF"="SYMBOL"))
    exp_data_input$ncbi_symbol<- exp_data_input$ID_REF
    colnames(exp_data_input)[colnames(exp_data_input)== "ENTREZID"]<- "ncbi_id"
  } else if(iden_type_list[iden_type_indx]=="ENTREID") {
    library(EnsDb.Hsapiens.v79)
    res_df <- ensembldb::select(EnsDb.Hsapiens.v79, keys= exp_data_input$ID_REF, keytype = "ENTREZID", columns = c("ENTREZID") )
    exp_data_input<- dplyr::left_join(exp_data_input, res_df, by= c("ID_REF"="ENTREZID"))
    exp_data_input$ncbi_id<- exp_data_input$ID_REF
    colnames(exp_data_input)[colnames(exp_data_input)== "SYMBOL"]<- "ncbi_symbol"
  } else{
    library(EnsDb.Hsapiens.v79)
    res_df <- ensembldb::select(EnsDb.Hsapiens.v79, keys= exp_data_input$ID_REF, keytype = iden_type_list[iden_type_indx], columns = c("SYMBOL","ENTREZID") )
    res_df <- res_df[!duplicated(res_df[,iden_type_list[iden_type_indx]]),]
    exp_data_input$ncbi_id = res_df$ENTREZID
    exp_data_input$ncbi_symbol = res_df$SYMBOL
  }
  #combining extracts IDs and Symbols to filtered data occurs in each statements of identifier conversion

  #Removal of ncbi_id with NA values as it represents no relevant gene
  exp_data_input <- dplyr::filter(exp_data_input, ncbi_symbol != "NA")
  exp_data_input <- dplyr::filter(exp_data_input, ncbi_id != "NA")

  # resolve duplicated genes if found
  if (sum(duplicated(exp_data_input$ncbi_symbol)) >0) {
    library(dplyr)
    exp_data_input <- exp_data_input %>% group_by(ncbi_symbol) %>% summarise_all(max)
  }

  #Check repeated individual samples in defined groups
  #test loop to extract individual sample names to be removed
  if (sum(sample_data_input$individuals == sample_data_input$accession)==0){

    samp_clean <- dplyr::group_by(sample_data_input, sample_data_input$groups) %>% summarize(which(duplicated(individuals)))

    dup_index= samp_clean$`which(duplicated(individuals))`

    dup_gp = samp_clean$'sample_data_input$groups'

    if(length(dup_index)!=0){
      dup_samp=c()
      for (x in 1:length(dup_index)){
        temp_sam_df= dplyr::filter(sample_data_input, sample_data_input$groups==dup_gp[x])
        dup_samp= c(dup_samp, temp_sam_df$individuals[dup_index[x]])
        sample_data_input= dplyr::filter(sample_data_input, !individuals %in% dup_samp)#remove duplicated sample from sample data
      }}
  }


  #update expression file and drop ID_REF column
  sample_data_input <- data.frame(sample_data_input)
  new_accn <- sample_data_input$accession
  exp_data_input <- dplyr::select(exp_data_input, c(all_of(new_accn), "ncbi_id", "ncbi_symbol"))

  #saving data for further filtration
  exp_data <<- exp_data_input
  sample_data <<- sample_data_input

  # Preparing data for analysis
  exp_data_input <- data.frame(exp_data_input)
  rownames(exp_data_input) <- exp_data_input$ncbi_symbol
  exp_data_input <<- exp_data_input[,sample_data$accession]

  rownames(sample_data_input) <- sample_data_input$accession
  sample_data_input <<- sample_data_input[, c("groups", "individuals")]

}

#'Function set number of Up- and down- regulated genes and build visual analysis: Pearson's correlation matrix for each group; Clustered heatmap and annotated to 20 pathways in both Up and Down regulated
#' @export

vis_function <- function(){
  #call library

  #adjust top DEG gene number for analysis

  print("Please select number of DEGs for analysis and building visulaizations step")

  deg_sel_menu <- c("30 up: 0 down", "25 up: 5 down", "20 up: 10 down", "15 up: 15 down",
                    "10 up: 20 down", "5 up: 25 down", "0 up: 30 down")

  #setting up regulated number; down regulated can be calculated as 40 - this value
  deg_up_number <- c(30, 25, 20, 15, 10, 5, 0)
  deg_index <- menu(deg_sel_menu, title= "DEG count to be analysed", graphics = TRUE)
  #set selected number of genes to be extracted
  deg_up_n <- as.numeric(deg_up_number[deg_index])
  deg_down_n <- 30 - as.numeric(deg_up_number[deg_index])
  #filter Analysis results for significant DEG
  filtered_df<-dplyr::filter(DEG_result_table,padj<adj_threshold)
  if (dim(filtered_df)[1] ==0){
    print("no result output to analysis")
  }
  #get up regulated gene names
  gene_id_up=c()
  for (i in 1:deg_up_n){
    gene_id_up= c(gene_id_up, filtered_df$ncbi_symbol[order(filtered_df$log2FoldChange, decreasing=TRUE)][i])
  }
  #get down regulated gene names
  gene_id_down=c()
  for (i in 1:deg_down_n){
    gene_id_down= c(gene_id_down, filtered_df$ncbi_symbol[order(filtered_df$log2FoldChange, decreasing=FALSE)][i])
  }
  #combine up- and down- regulated gene_ID
  gene_id <- c(gene_id_up, gene_id_down)

  #Filtered expression data frame for groups
  merged_exp_deg<- exp_data_input[gene_id,total_sample_names]
  group_1_exp_deg <-exp_data_input[gene_id,group_1_samples]
  group_2_exp_deg <-exp_data_input[gene_id,group_1_samples]


  library(ggcorrplot)

  print("If a standard deviation is zero appears: Gene with zero R value include same repeated value most probably zero")

  library(pheatmap)
  library (matrixStats)

  #working on merged sata of both  DEGs files of normal and primary cancer
  merged_exp_deg_log2<- log2(merged_exp_deg+1) #get the log data
  merged_exp_deg__meansubtract<- merged_exp_deg_log2-rowMeans((merged_exp_deg_log2))#get the mean of data
  merged_exp_deg_zscores<- merged_exp_deg__meansubtract/rowSds(as.matrix(merged_exp_deg_log2))#get the z-score

  #clustering
  #eucledian_distance
  #compelete_linkage
  merged_exp_deg__meansubtract%>% dist() #eucledian_distance
  merged_exp_deg_rowclust<- merged_exp_deg__meansubtract%>% dist()%>%hclust()#clustering, compelete_linkage

  ##correlation
  merged_cor_metrix<- cor(t(merged_exp_deg__meansubtract))#correlation matrix for combined data
  merged_dist <- as.dist(1-merged_cor_metrix) #eucledian_distance
  merged_rowclust_cor= hclust(merged_dist)#clustering, compelete_linkage

  #annotation
  annotation_col= sample_data_input[1] #make annotation for columns
  data<- filtered_df[,c("log2FoldChange", "ncbi_symbol", "ncbi_id")]
  rownames(data) <- filtered_df[,c("ncbi_symbol")]#for rows, get the data from filtered df
  data3<- data[rownames(data)%in% gene_id,]#filter the file with our chosen genes
  data3$Gene_class=data3[,c('log2FoldChange')] #duplicate the foldchange column to new column
  data3$Gene_class [data3$Gene_class  > 0] <- "UP_DEGs" #replace the new column with up genes
  data3$Gene_class [data3$Gene_class  < 0] <- "Down_DEGs" ##replace the new column with down genes
  annotation_row=data3["Gene_class"] #select the gene class column for annotation


  setwd(data_dir)
  write.csv(cor(na.omit(t(group_1_exp_deg))), "Group_1_correlation.csv")
  write.csv(cor(na.omit(t(group_2_exp_deg))), "Group_2_correlation.csv")

  setwd(vis_dir)

  png("Group_1_corr.png", res=300, width=7, height=4.5, unit="in")
  plot(ggcorrplot::ggcorrplot(cor(t(group_1_exp_deg)),type= "lower",
                              tl.cex= 7, tl.srt=90,lab = TRUE ,lab_size = 1,
                              title= group_1_name , show.diag = TRUE,
                              legend.title = "Pearson's \n corr"))

  dev.off()

  png("Group_2_corr.png", res=300, width=7, height=4.5, unit="in")
  plot(ggcorrplot::ggcorrplot(cor(na.omit(t(group_2_exp_deg))),type= "lower",
                              tl.cex= 7, tl.srt=90, lab = TRUE ,lab_size = 1,
                              title= group_2_name, show.diag = TRUE,
                              legend.title = "Pearson's \n corr"))
  dev.off()

  png("Clustered_heatmap.png", res=300, width=7, height=4.5, unit="in")
  pheatmap(merged_exp_deg__meansubtract, clustering_distance_rows = "correlation", cluster_cols = T, annotation_col = annotation_col,
           annotation_row = annotation_row, main = "Clustered heatmap",  fontsize=5)


  dev.off()

  print("Pathway extraction starting...")
  ###### KEGG pathways
  library(pathview)
  library(gage)
  library(gageData)

  kegg_foldchanges <- filtered_df$log2FoldChange
  names(kegg_foldchanges) <- filtered_df$ncbi_id

  data("kegg.sets.hs")
  #signaling metabolic pathways
  data("sigmet.idx.hs")

  kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

  keggres <- gage(kegg_foldchanges, gsets= kegg.sets.hs, same.dir=TRUE)

  library(dplyr)
  #creating up regulated data frame
  Keggrespathways <- data.frame(id= rownames(keggres$greater), keggres$greater ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(row_number() <=20) %>%
    .$id %>%
    as.character()

  keggresid <- substr(Keggrespathways, start=1, stop=8)

  #change path
  setwd(new_dir_up)

  # Plotting up regulated pathways
  tmp_plot= sapply(keggresid, function(pid) pathview(gene.data= kegg_foldchanges, pathway.id= pid, species="hsa"))

  # Creating down regulated genes data frame
  Keggrespathways = data.frame(id= rownames(keggres$less), keggres$less ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(row_number() <=20) %>%
    .$id %>%
    as.character()

  keggresid =substr(Keggrespathways, start=1, stop=8)
  #Change path

  setwd(new_dir_down)

  #Plotting down regulated pathways
  tmp_plot= sapply(keggresid, function(pid) pathview(gene.data= kegg_foldchanges, pathway.id= pid, species="hsa"))
}

#'Function define two groups, perform analysis using DESeq2 library and visulaize data with set padj and log fold change as a Volcano plot
#' @export

deseq_analysis <- function(){
  #call library

  #Define groups for analysis
  defined_groups <- unique(sample_data_input$groups)
  #group 1
  group_1_index<-menu(defined_groups, title= "Define first group for analysis", graphics=TRUE)
  group_1_name<<- defined_groups[group_1_index]
  #group 2
  group_2_index<-menu(defined_groups, title= "Define second group for analysis", graphics=TRUE)
  group_2_name<<- defined_groups[group_2_index]
  # Checkpoint for groups
  if ((group_1_name == group_2_name)) {
    print("Please define two different groups")
    break
  } else if(length(defined_groups) <2){
    print("No enough defined groups; Please re-read the sample data and select column of group definition")
    break
  } else if(length(defined_groups) ==2){
    group_1_samples<<- rownames(sample_data_input[sample_data_input$groups ==defined_groups[group_1_index],])
    group_2_samples<<- rownames(sample_data_input[sample_data_input$groups ==defined_groups[group_2_index],])
    total_sample_names <<- c(group_1_samples, group_2_samples)
    analysis_pass_genes <- rownames(exp_data_input)
  } else{ #for more than two groups; genes should be re-filtered for low expression
    group_1_samples<<- rownames(sample_data_input[sample_data_input$groups ==defined_groups[group_1_index],])
    group_2_samples<<- rownames(sample_data_input[sample_data_input$groups ==defined_groups[group_2_index],])
    total_sample_names <<- c(group_1_samples, group_2_samples)
    #filtering genes of low expression in these two groups
    #This step was done in cleaning filtration step for all samples
    analysis_pass_genes <- rownames(dplyr::filter(exp_data_input, rowSums(exp_data_input[,total_sample_names]) > ((length(total_sample_names))*2)))
  } #end of >2 if statement


  #change working directory
  setwd(data_dir)
  write.csv (exp_data, "Exp_analysis_input.csv")
  write.csv (sample_data, "Sample_data_input.csv")
  write.csv (abs(round(exp_data_input[analysis_pass_genes,total_sample_names], 0)),"groups_filt_expression.csv")
  write.csv(sample_data_input[total_sample_names,], "groups_filt_sample_data.csv")

  library(DESeq2)

  #Ordering columns and rows in the two files
  exp_data_input <- exp_data_input[,rownames(sample_data_input)]

  deseq2Data <- DESeqDataSetFromMatrix(countData = abs(round(exp_data_input[analysis_pass_genes,total_sample_names], 0)),
                                       colData = sample_data_input[total_sample_names,],
                                       design = ~ groups)


  print("Analysis started")
  #Differential Expression Analysis
  deseq2Data <- DESeq(deseq2Data)

  # Extract differential expression results
  # For "tissueType" [ primary vs normal] comparison
  deseq2Results <- results(deseq2Data,
                           contrast=c("groups",
                                      group_1_name,
                                      group_2_name))


  #####DESeq2 results filtration and cleaning
  res_df<<-as.data.frame(deseq2Results) #reading results into dataframe

  print("Analysis done!")
  #Merge entrez id to res_df
  res_df$ncbi_symbol<- rownames(res_df)
  ext_list <- rownames(res_df)
  exp_mod<- exp_data[exp_data$ncbi_symbol %in% ext_list,]
  res_df<- dplyr::left_join(res_df, exp_mod[,c("ncbi_id", "ncbi_symbol")], by= "ncbi_symbol")

  #set adj. p. value threshold
  adj_p_val_threshold <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)

  print(paste0("Min adjusted pvalue calculated:     ",  min(na.omit(res_df$padj))))

  print(paste0("Max adjusted pvalue calculated:     ",  max(na.omit(res_df$padj))))

  print("Please select adj. p. value threshold (FDR)")

  adj_threshold_indx <- menu(as.character(adj_p_val_threshold), title= "adj.p.value", graphics= TRUE)

  adj_threshold<<- adj_p_val_threshold[adj_threshold_indx]


  print("Please select log fold change cutoff value")

  lfc_threshold_list <- c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, "default = 0.8")

  print(paste0("Min log2 fold change calculated:     ", min(na.omit(res_df$log2FoldChange))))
  print(paste("Max log2 fold change calculated:     ", max(na.omit(res_df$log2FoldChange))))


  lfc_threshold_indx <- menu(as.character(lfc_threshold_list), title= "log fold change threshold", graphics= TRUE)

  if (lfc_threshold_indx == 10){
    lfc_threshold <<- 0.8
  } else{
    lfc_threshold <<- as.numeric(lfc_threshold_list[lfc_threshold_indx])
  }

  ##### visualization

  library(ggrepel)
  library(ggplot2)

  #Volcano pot
  res_df_volcano<- na.omit(res_df)
  res_df_volcano$diffexpression <- "NO"
  #labeling up regulated and down regulated genes
  res_df_volcano$diffexpression[res_df_volcano$log2FoldChange >= lfc_threshold
                                & res_df_volcano$padj<=adj_threshold ] <-"UP regulated"
  res_df_volcano$diffexpression[res_df_volcano$log2FoldChange< lfc_threshold
                                & res_df_volcano$padj<adj_threshold ] <-"Down regulated"

  # Labelling top 10 most significant genes
  res_df_volcano$delabel <- NA
  gen_lab_up_threshold<- (arrange(res_df_volcano, pvalue))$pvalue[10]
  gen_lab_up <- (arrange(res_df_volcano, pvalue))$ncbi_symbol[1:10]

  res_df_volcano$delabel[res_df_volcano$pvalue<= gen_lab_up_threshold]<- (res_df_volcano$ncbi_symbol[res_df_volcano$pvalue<=gen_lab_up_threshold])
  res_df_volcano<<- res_df_volcano

  #exporting plot
  setwd(data_dir)
  write.csv(res_df_volcano[,1: ncol(res_df_volcano)-1], "DEG_analysis_results.csv")

  setwd(vis_dir)

  pdf("DESeq2_Volcano_plot.pdf")
  plot(ggplot(data= res_df_volcano, aes(x=log2FoldChange, y= -log10(pvalue), col= diffexpression, label=delabel))+
         geom_point()+
         theme_minimal()+
         geom_text_repel()+
         scale_color_manual(values= c("blue", "black", "red"))+
         geom_vline(xintercept = c(-lfc_threshold, lfc_threshold))+
         geom_hline(yintercept = -log10(adj_threshold))+
         theme(text=element_text(size=11))+
         ggtitle("Differential expression in genes with top 10 significant genes labelled"))

  dev.off()

  #exporting data frames
  expression_data_analyzed <<- exp_data_input[analysis_pass_genes,total_sample_names]
  sample_data_analyzed <<- sample_data_input[total_sample_names,]
  DEG_result_table<<- res_df_volcano[,1: ncol(res_df_volcano)-1]

}

#'Function check for library dependencies (first time) then reads user files as expression data and sample phenotype data as tab delimited files and extract needed data for further analysis; format inputs for data_filt(), deseq_functiona nd vis-function
#' @export

read_external_file <- function (){

  #check libraries for first time- this section will be skipped automatically after first use
  if (!"BiocManager" %in%installed.packages() &&
      BiocManager::version() != 3.15){
    install.packages("BiocManager")
  }

  if (!"remotes" %in% installed.packages()){
    install.packages("remotes")
  }


  ### checking for packages and install dependencies
  pack_list<- c("pheatmap","ggcorrplot", "ggrepel", "ggplot2", "matrixStats", "dplyr", "stringr",
                "R.utils", "tibble" )
  pack_ver <- c("1.0.12", "0.1.4","0.9.1", "3.3.6", "0.62.0", "1.0.10", "1.4.1", "2.12.0", "3.1.8")
  required_packs <- c()
  required_p_ver <- c()

  for (i in 1: length(pack_list)) {
    if (!pack_list[i] %in% installed.packages()){
      required_packs <- c(required_packs, pack_list[i])
      required_p_ver <- c(required_p_ver, pack_ver[i])
    }
  }

  if (length(required_packs>0)){
    for (i in 1: length(required_packs)){
      remotes::install_version(required_packs[i], version = required_p_ver, repos = "http://cran.us.r-project.org")
    }
  }

  # checking for biconductor dependencies

  b_pack_list<- c("GEOquery","gageData", "gage", "pathview", "SummarizedExperiment", "MatrixGenerics",
                  "EnsDb.Hsapiens.v79","ensembldb", "AnnotationFilter", "GenomicFeatures", "AnnotationDbi",
                  "Biobase", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "BiocGenerics",
                  "ensembldb", "illuminaHumanv4.db", "org.Hs.eg.db", "DESeq2" )
  b_pack_ver <- c("2.64.2","2.34.0", "2.46.1" , "1.36.1", "1.26.1", "1.8.1", "2.99.0", "2.20.2",
                  "1.20.0", "1.48.4", "1.58.0", "2.56.0", "1.48.0", "1.32.4", "2.30.1", "0.34.0", "0.42.0",
                  "2.20.2", "1.26.0" , "3.15.0", "1.36.0" )
  b_required_packs <- c()
  b_required_p_ver <- c()

  for (i in 1: length(b_pack_list)) {
    if (!b_pack_list[i] %in% installed.packages()){
      b_required_packs <- c(b_required_packs, b_pack_list[i])
      b_required_p_ver <- c(b_required_p_ver, b_pack_ver[i])
    }
  }

  if (length(b_required_packs) >0){
    for (i in 1: length(b_required_packs)){
      BiocManager::install(b_required_packs[i], version = b_required_p_ver[i])
    }
  }

  library(dplyr)
  #creating job directory
  curr_path<<- getwd()
  job_dir<<- paste0(curr_path, "/exp2vis_analysis")

  if (dir.exists(job_dir)){
    print("exp2vis_analysis folder exists! please move or rename folder then continue")
    dir_check<- menu(c("Folder moved, continue", "cancel"), title = "Please rename/ move exp2vis_analysis folder", graphics = TRUE)
  } else {
    dir_check = 0
  }

  if (dir_check ==2){
    break
  }

  dir.create(job_dir)

  #creating data files directory
  data_dir <<- paste0(job_dir, "/data_files")
  dir.create(data_dir)

  #creating visualization directory
  vis_dir<<- paste0(job_dir, "/visualization")
  dir.create(vis_dir)

  #creating up regulated pathway directory
  new_dir_up<<- paste0(job_dir, "/upregulated_DEG_pathways")
  dir.create(new_dir_up)

  #creating up regulated pathway directory
  new_dir_down<<- paste0(job_dir, "/down-regulated_DEG_pathways")
  dir.create(new_dir_down)


  print("Please select sample data file (Phenotype) as tab delimited text" )
  #get file name
  s_file_index <- menu(list.files(getwd()), title= "Select sample data file", graphics=TRUE)

  sample_file <- list.files(getwd())[s_file_index]
  #read file
  sample_data <- read.delim(file =sample_file, sep= "\t")

  #Extract sample data
  View(sample_data[1:5,])
  print("Please select column with sample identifier")
  #index for accessions
  accn_col_def <- menu(colnames(sample_data), title= "Sample identifier", graphics= TRUE)

  print("Please select column number to define groups for analysis; if not found please cancel")
  #index for group definition
  group_col_def <-menu(colnames(sample_data), title= "Select column to define groups", graphics= TRUE)

  print("Please select column number to define individual samples; e.g : donors, individual id")
  print("If no column represents individuals please selection column for GEO accession")

  #index for individual identifier
  indiv_col_def <-menu(colnames(sample_data), title= "Select column number to define groups for analysis; if not found please cancel", graphics= TRUE)

  #Create sample data file
  sample_data_input <- sample_data[, c(accn_col_def,group_col_def, indiv_col_def )]
  colnames(sample_data_input) <- c("accession", "groups", "individuals")

  # Expression file
  print("Please select expression data file as tab delimited text; make sure to include gene identifier under column name ID_REF" )

  #get file name
  ex_file_index <- menu(list.files(getwd()), title= "Select expression data file", graphics=TRUE)

  exp_file <- list.files(getwd())[ex_file_index]
  #read file
  exp_data <- read.delim(file =exp_file, sep= "\t")

  exp_data_col <- c("ID_REF", sample_data_input$accession)
  samples_iden_cols <- exp_data_col[-which(exp_data_col =="ID_REF")]
  setequal(samples_iden_cols, sample_data_input$accession)

  if (setequal(colnames(exp_data), exp_data_col) ==FALSE) {
    print("Please remove columns other than samples and ID_REF from expression file")
  } else if(setequal(samples_iden_cols, sample_data_input$accession)==FALSE){
    print(" Samples identifier in two files does not match; only found in expression file will be kept")
    #code should be added to filter sample data from unused samples in expression file
  }
  exp_data_input <<- exp_data

  sample_data_input <<- sample_data_input

  data_filt()

  #analysis checkpoint

  print("please select analysis pipeline")
  package_sel<- c("DESeq2 package: NB: package rounds up expression values to zero decimels ", "Limma")

  analysis_index <- menu (package_sel, title= "Select analysis pipeline", graphics =TRUE)

  if (analysis_index ==1){
    deseq_analysis()
  } else {
    limma_analysis()
  }

  print("Analysis is all done! building Visuals step starting")

  vis_function()


}


#'Function check for library dependencies (first time) then reads GEO GSE ID; specifies platform -GPL- extracts expression and pheotype data then format data inputs for data_filt()m deseq_function(),and vis_function
#' @export
#' @param gse_id character of GEO GSE ID

read_gse <- function(gse_id){

  #check libraries for first time- this section will be skipped automatically after first use
  if (!"BiocManager" %in%installed.packages() &&
      BiocManager::version() != 3.15){
    install.packages("BiocManager")
  }

  if (!"remotes" %in% installed.packages()){
    install.packages("remotes")
  }


  ### checking for packages and install dependencies
  pack_list<- c("pheatmap","ggcorrplot", "ggrepel", "ggplot2", "matrixStats", "dplyr", "stringr",
                "R.utils", "tibble" )
  pack_ver <- c("1.0.12", "0.1.4","0.9.1", "3.3.6", "0.62.0", "1.0.10", "1.4.1", "2.12.0", "3.1.8")
  required_packs <- c()
  required_p_ver <- c()

  for (i in 1: length(pack_list)) {
    if (!pack_list[i] %in% installed.packages()){
      required_packs <- c(required_packs, pack_list[i])
      required_p_ver <- c(required_p_ver, pack_ver[i])
    }
  }

  if (length(required_packs>0)){
    for (i in 1: length(required_packs)){
      remotes::install_version(required_packs[i], version = required_p_ver, repos = "http://cran.us.r-project.org")
    }
  }

  # checking for biconductor dependencies

  b_pack_list<- c("GEOquery","gageData", "gage", "pathview", "SummarizedExperiment", "MatrixGenerics",
                  "EnsDb.Hsapiens.v79","ensembldb", "AnnotationFilter", "GenomicFeatures", "AnnotationDbi",
                  "Biobase", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "BiocGenerics",
                  "ensembldb", "illuminaHumanv4.db", "org.Hs.eg.db", "DESeq2" )
  b_pack_ver <- c("2.64.2","2.34.0", "2.46.1" , "1.36.1", "1.26.1", "1.8.1", "2.99.0", "2.20.2",
                  "1.20.0", "1.48.4", "1.58.0", "2.56.0", "1.48.0", "1.32.4", "2.30.1", "0.34.0", "0.42.0",
                  "2.20.2", "1.26.0" , "3.15.0", "1.36.0" )
  b_required_packs <- c()
  b_required_p_ver <- c()

  for (i in 1: length(b_pack_list)) {
    if (!b_pack_list[i] %in% installed.packages()){
      b_required_packs <- c(b_required_packs, b_pack_list[i])
      b_required_p_ver <- c(b_required_p_ver, b_pack_ver[i])
    }
  }

  if (length(b_required_packs) >0){
    for (i in 1: length(b_required_packs)){
      BiocManager::install(b_required_packs[i], version = b_required_p_ver[i])
    }
  }



  #creating job directory
  curr_path<<- getwd()
  job_dir<<- paste0(curr_path, "/exp2vis_analysis")

  if (dir.exists(job_dir)){
    print("exp2vis_analysis folder exists! please move or rename folder then continue")
    dir_check<- menu(c("Folder moved, continue", "cancel"), title = "Please rename/ move exp2vis_analysis folder", graphics = TRUE)
  } else {
    dir_check = 0
  }

  if (dir_check ==2){
    break
  }

  dir.create(job_dir)

  #creating data files directory
  data_dir <<- paste0(job_dir, "/data_files")
  dir.create(data_dir)

  #creating visualization directory
  vis_dir<<- paste0(job_dir, "/visualization")
  dir.create(vis_dir)

  #creating up regulated pathway directory
  new_dir_up<<- paste0(job_dir, "/upregulated_DEG_pathways")
  dir.create(new_dir_up)

  #creating up regulated pathway directory
  new_dir_down<<- paste0(job_dir, "/down-regulated_DEG_pathways")
  dir.create(new_dir_down)


  library(GEOquery)

  #set GSE ID from input parameter AnnotGPL=TRUE
  #Read GSE data matrix
  gse <- getGEO(gse_id,GSEMatrix=TRUE)

  if (length(gse)>1){
    platforms_sel=c()
    for (i in 1: length(gse)){
      platforms_sel= c(platforms_sel, gse[[i]]@annotation)
    }
    print("Please selcect platform")

    gse_idx<- menu(platforms_sel, title= "Select platform GPLxxx", graphics = TRUE)
  } else{
    gse_idx = 1}

  #expression data extraction
  exp_df<- data.frame(gse[[gse_idx]]@assayData$exprs)

  if(dim(exp_df)[1] ==0){
    print("Expression data empty!")

    break
  }

  #some expression data is acquired with identifiers rownames while others are just numbered
  #Converting rownames to numbers to avoid confusion
  rownames(exp_df) <- seq(1, dim(exp_df)[1])

  #read feature data
  gse_fd<- gse[[gse_idx]]@featureData@data
  #Converting rownames to numbers to avoid confusion
  rownames(gse_fd) <- seq(1, dim(gse_fd)[1])

  #define first 10 rowns in phenotype data
  gse_fd_head<- head (gse_fd, 10)

  #Acquire ncbi symbol from data
  print("please identify column including NCBI symbol identifier e.g. DDR1")
  View(gse_fd_head)
  gse_iden_index <- menu(colnames(gse_fd_head), title = "NCBI symbol column", graphics = TRUE)

  #extracting column to vector
  ID_REF <- gse_fd[,gse_iden_index]
  #cleaning in case of more than one record in the same value e.g DDR1 /// MIR4640
  ID_REF <- sub(" .*","" ,ID_REF)

  #Extract gene identifier
  exp_df$ID_REF <- ID_REF

  #Extract sample data

  gse_pd_df <- data.frame(gse[[gse_idx]]@phenoData@data [1:10,])
  #change rownames
  rownames(gse_pd_df) <- seq(1, dim(gse_pd_df)[1])

  View(gse_pd_df)
  #Define sample GEO accession
  print("Please define group with sample GEO accession identifier GSMxxx")

  accn_col_def<-menu(colnames(gse_pd_df), title= "Select column to define GEO accession", graphics=TRUE)


  #Define groups columns
  print("Please select column number to define groups for analysis; Please wait till dataframe uploads. if not found please cancel")

  group_col_def<-menu(colnames(gse_pd_df), title= "Select column to define groups", graphics=TRUE)

  if (group_col_def==0 ){
    print( "Sample data can not to retrieved; please reload the function and select another format")
    break
  }

  #Define individuals column
  print("Please select column number to define individual samples; e.g : donors, individual id")
  print("If no column represents individuals please selection column for GEO accession")

  indiv_col_def<-menu(colnames(gse[[gse_idx]]@phenoData@data), title= "Select column to define individuals", graphics=TRUE)

  #Create sample dataframe
  sample_data_input <- data.frame (matrix(ncol = 3, nrow = length(gse[[gse_idx]]@phenoData@data$geo_accession)))
  colnames(sample_data_input) <- c("accession", "groups", "individuals")
  sample_data_input$accession =gse[[gse_idx]]@phenoData@data[, accn_col_def]
  sample_data_input$groups =gse[[gse_idx]]@phenoData@data[,group_col_def]
  sample_data_input$individuals= gse[[gse_idx]]@phenoData@data[,indiv_col_def]

  # exp_data_col <<- c("ID_REF", sample_data_input$accession)
  # samples_iden_cols <<- exp_data_col[-which(exp_data_col =="ID_REF")]

  sample_data_input <<- sample_data_input
  exp_data_input <<- exp_df



  data_filt()

  #analysis checkpoint

  print("please select analysis pipeline")
  package_sel<- c("DESeq2 package: NB: package rounds up expression values to zero decimels ", "Limma")

  analysis_index <- menu (package_sel, title= "Select analysis pipeline", graphics =TRUE)

  if (analysis_index ==1){
    deseq_analysis()
  } else {
    limma_analysis()
  }

  print("Analysis is all done! building Visuals step starting")

  vis_function()


}


#'Function check for library dependencies (first time) then reads GEO GSE series matrix zipped file; extracts expression and pheotype data then format data inputs for data_filt()m deseq_function(),and vis_function
#' @export

read_series_matrix<- function (){

  #check libraries for first time- this section will be skipped automatically after first use
  if (!"BiocManager" %in%installed.packages() &&
      BiocManager::version() != 3.15){
    install.packages("BiocManager")
  }

  if (!"remotes" %in% installed.packages()){
    install.packages("remotes")
  }


  ### checking for packages and install dependencies
  pack_list<- c("pheatmap","ggcorrplot", "ggrepel", "ggplot2", "matrixStats", "dplyr", "stringr",
                "R.utils", "tibble" )
  pack_ver <- c("1.0.12", "0.1.4","0.9.1", "3.3.6", "0.62.0", "1.0.10", "1.4.1", "2.12.0", "3.1.8")
  required_packs <- c()
  required_p_ver <- c()

  for (i in 1: length(pack_list)) {
    if (!pack_list[i] %in% installed.packages()){
      required_packs <- c(required_packs, pack_list[i])
      required_p_ver <- c(required_p_ver, pack_ver[i])
    }
  }

  if (length(required_packs>0)){
    for (i in 1: length(required_packs)){
      remotes::install_version(required_packs[i], version = required_p_ver, repos = "http://cran.us.r-project.org")
    }
  }

  # checking for biconductor dependencies

  b_pack_list<- c("GEOquery","gageData", "gage", "pathview", "SummarizedExperiment", "MatrixGenerics",
                  "EnsDb.Hsapiens.v79","ensembldb", "AnnotationFilter", "GenomicFeatures", "AnnotationDbi",
                  "Biobase", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "BiocGenerics",
                  "ensembldb", "illuminaHumanv4.db", "org.Hs.eg.db", "DESeq2" )
  b_pack_ver <- c("2.64.2","2.34.0", "2.46.1" , "1.36.1", "1.26.1", "1.8.1", "2.99.0", "2.20.2",
                  "1.20.0", "1.48.4", "1.58.0", "2.56.0", "1.48.0", "1.32.4", "2.30.1", "0.34.0", "0.42.0",
                  "2.20.2", "1.26.0" , "3.15.0", "1.36.0" )
  b_required_packs <- c()
  b_required_p_ver <- c()

  for (i in 1: length(b_pack_list)) {
    if (!b_pack_list[i] %in% installed.packages()){
      b_required_packs <- c(b_required_packs, b_pack_list[i])
      b_required_p_ver <- c(b_required_p_ver, b_pack_ver[i])
    }
  }

  if (length(b_required_packs) >0){
    for (i in 1: length(b_required_packs)){
      BiocManager::install(b_required_packs[i], version = b_required_p_ver[i])
    }
  }


  #creating job directory
  curr_path<<- getwd()
  job_dir<<- paste0(curr_path, "/exp2vis_analysis")

  if (dir.exists(job_dir)){
    print("exp2vis_analysis folder exists! please move or rename folder then continue")
    dir_check<- menu(c("Folder moved, continue", "cancel"), title = "Please rename/ move exp2vis_analysis folder", graphics = TRUE)
  } else {
    dir_check = 0
  }

  if (dir_check ==2){
    break
  }

  dir.create(job_dir)

  #creating data files directory
  data_dir <<- paste0(job_dir, "/data_files")
  dir.create(data_dir)

  #creating visualization directory
  vis_dir<<- paste0(job_dir, "/visualization")
  dir.create(vis_dir)

  #creating up regulated pathway directory
  new_dir_up<<- paste0(job_dir, "/upregulated_DEG_pathways")
  dir.create(new_dir_up)

  #creating up regulated pathway directory
  new_dir_down<<- paste0(job_dir, "/down-regulated_DEG_pathways")
  dir.create(new_dir_down)

  #user input for downloaded series matrix file
  print( "Please make sure that series matrix zipped file is copied in the working directory")

  file_index<-menu(list.files(getwd()), title= "Select Zipped matrix file", graphics=TRUE)

  zip_file_name= list.files(getwd())[file_index]
  #zip_file_name=
  #unpack zipped file
  R.utils::gunzip(zip_file_name)

  print("File unpacked successfully ")

  txt_file_name = stringr::str_replace(zip_file_name, ".gz", "") #create text file name
  #Expression matrix extraction
  mat_file <- read.delim(file =txt_file_name, sep= "\t", skip=which(readLines(txt_file_name)=="!series_matrix_table_begin"),comment.char = "!")

  if(dim(mat_file)[1] ==0){
    print("No expression data retreived; please reload function and try another format")
    break
  }

  #Sample data extraction
  sub_info_line = which(grepl("!Series_relation", readLines(txt_file_name)))
  first_line <- read.table(file= txt_file_name, sep="\t", fill=TRUE, skip=sub_info_line+1)
  sample_data_ext<- t(first_line)
  colnames(sample_data_ext)<- seq(1, ncol(sample_data_ext))
  sample_data_ext<- sample_data_ext[-1,]

  View(head(sample_data_ext))

  print("Please select column number for sample accession (GEO accession")
  accn_col_def<-menu(colnames(sample_data_ext), title= "Select sample accession", graphics=TRUE)


  print("Please select column number to define groups for analysis; if not found please cancel")
  group_col_def<-menu(colnames(sample_data_ext), title= "Select column to define groups", graphics=TRUE)

  if (group_col_def==0 ){
    cat( "Sample data can not to retrieved; please reload the function and select another format")
    break
  }




  print("Please select column number to define individual samples; e.g : donors, individual id")
  print("If no column represents individuals please selection column for GEO accession")
  indiv_col_def<-menu(colnames(sample_data_ext), title= "Select column to define individuals", graphics=TRUE)

  sample_data_input <- sample_data_ext[,c(accn_col_def, group_col_def, indiv_col_def)]
  colnames(sample_data_input) <- c("accession", "groups", "individuals")

  sample_data_input <<- data.frame(sample_data_input)
  exp_data_input <<- mat_file


  data_filt()

  #analysis checkpoint

  print("please select analysis pipeline")
  package_sel<- c("DESeq2 package: NB: package rounds up expression values to zero decimels ", "Limma")

  analysis_index <- menu (package_sel, title= "Select analysis pipeline", graphics =TRUE)

  if (analysis_index ==1){
    deseq_analysis()
  } else {
    limma_analysis()
  }

  print("Analysis is all done! building Visuals step starting")

  vis_function()

}


