library(shiny)
library(magrittr)
library(ggplot2)
library(data.table)
library(ggpubr)

# Server Logic
shinyServer(function(input, output) {
    
    # ----
    # Load CCLE Cluster Data
    # ----
    cl_path <- "CCLE/Clusters/"
    cl2_files <- list.files(cl_path)
    cluster_filebase <- cl2_files[grepl("csv", cl2_files) & grepl("gmm", cl2_files)] %>% tools::file_path_sans_ext()
    cluster_filepath <- paste0(cl_path, "/", cluster_filebase, ".csv")
    
    clusters_raw <- list()
    for(i in 1:length(cluster_filepath)){
        clusters_raw[[i]] <- data.table::fread(cluster_filepath[i], data.table = FALSE)
        names(clusters_raw)[i] <- cluster_filebase[i]
    }
    
    # Annotate with Tissue Type
    ccle_annot_path <- "CCLE/sample_info.csv"
    ccle_annot <- data.table::fread(ccle_annot_path, data.table = FALSE)
    
    clusters_raw <- lapply(clusters_raw, FUN = function(x){
        colnames(x) <- c("Cell_Line", "Group")
        x$Tissue_Type <- ccle_annot$lineage[match(x$Cell_Line, ccle_annot$DepMap_ID)]
        x
    })
    
    cluster_ids <- c(
        "gmm_exp_2_pc_7_clusters.csv",
        "gmm_exp_55_pc_6_clusters.csv",
        "gmm_exp_55_pc_9_clusters.csv",
        "gmm_exp_snv_9_clusters.csv",
        "gmm_exp_snv_mutsig_8_clusters.csv",
        "gmm_mutsig_3_pc_9_clusters.csv",
        "gmm_snv_2_pc_8_clusters.csv"
    )
    cluster_descriptions <- c(
        "7 GMM clusters based on top 2 principal components from expression data",
        "6 GMM clusters based on PC's that explain >=70% of variance in the expression data",
        "9 GMM clusters based on PC's that explain >=70% of variance in the expression data",
        "9 GMM clusters based on PC's that explain >=70% of variance in the expression and SNV data",
        "8 GMM clusters based on PC's that explain >=70% of variance in the expression, SNV, and mutational signature data",
        "9 GMM clusters based on top 3 PC's from mutational signature data",
        "8 GMM clusters based on top 2 PC's from SNV data"
    )
    
    output$this_clust_desc <- renderText(paste0(
        "Model Description: ",
        cluster_descriptions[input$clustid]
    ))
    
    
    # ----
    # Plot Agreement Tissue Type
    # ----
    clust_tx_gg_list <- list()
    
    for(i in 1:length(clusters_raw)){
        #count of samples per tissue type for this cluster
        tissue_types <- clusters_raw[[i]]$Tissue_Type %>% unique
        n_per_tissue <- numeric()
        for(j in 1:length(tissue_types)){
            n_per_tissue[j] <- sum(clusters_raw[[i]]$Tissue_Type == tissue_types[j])
            names(n_per_tissue)[j] <- tissue_types[j]
        }
        
        # count of samples per group
        groups <- clusters_raw[[i]]$Group %>% unique %>% sort
        n_per_group_tissue <- data.frame(
            Tissue_Type = character(), Group = character(), FractionOfTissue = numeric(), Count = numeric())
        this_row <- data.frame(Tissue_Type = 0, Group = 0, FractionOfTissue = 0, Count = 0)
        for(k in 1:length(groups)){
            for(j in 1:length(tissue_types)){
                this_row$Tissue_Type <- tissue_types[j]
                this_row$Group <- groups[k]
                this_row$Count <- clusters_raw[[i]][
                    clusters_raw[[i]]$Group == groups[k] & clusters_raw[[i]]$Tissue_Type == tissue_types[j],] %>% nrow
                this_row$FractionOfTissue <- this_row$Count / n_per_tissue[names(n_per_tissue) == tissue_types[j]]
                
                n_per_group_tissue %<>% rbind(this_row) 
            }
        }
        
        # plot
        clust_tx_gg_list[[i]] <- ggplot(n_per_group_tissue) +
            geom_bar(aes(x = Tissue_Type, fill = Tissue_Type, y = FractionOfTissue), stat = "identity") +
            facet_grid(.~Group) +
            ggtitle(names(clusters_raw)[i]) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(x = Tissue_Type, y = FractionOfTissue, label = Count))
    }
    
    # ----
    # Load TCGA Data
    # ----
    tcga_cl_path <- "TCGA/Clusters/"
    tcga_cl_files <- list.files(tcga_cl_path)
    tcga_cl_files <- tcga_cl_files[grep(".csv", tcga_cl_files)]
    tcga_cl_names <- tcga_cl_files %>% tools::file_path_sans_ext()
    tcga_cl_list <- list()
    
    for(i in 1:length(tcga_cl_files)){
        tcga_cl_list[[i]] <- data.table::fread(paste0(tcga_cl_path, "/", tcga_cl_files[i]), data.table = FALSE)
    }
    
    tcga_annot_path <- "TCGA/tcga_sample_lists/"
    tcga_annot_files <- list.files(tcga_annot_path)
    tcga_annot_names <- tcga_annot_files %>% gsub("_samples\\.txt", "", .)

    tcga_annot_list <- list()
    
    for(i in 1:length(tcga_annot_files)){
        tcga_annot_list[[i]] <- readLines(paste0(tcga_annot_path, "/", tcga_annot_files[i]))
        names(tcga_annot_list)[i] <- tcga_annot_names[i]
    }
    
    tcga_annot_list <- lapply(tcga_annot_list, as.data.frame)
    
    tcga_annot_df <- dplyr::bind_rows(tcga_annot_list, .id = "Tissue_Type")
    colnames(tcga_annot_df)[2] <- "Tumor_ID"

    # add tissue data to clusters
    tcga_cl_list <- lapply(tcga_cl_list, FUN = function(x){
        x$Tissue_Type <- tcga_annot_df$Tissue_Type[match(x$Sample, tcga_annot_df$Tumor_ID)]
        x
    })
    
    names(tcga_cl_list) <- tcga_cl_names
    
    # ----
    # Plot TCGA clustering by tissue type
    # ----
    tcga_clust_tx_gg_list <- list()
    for(i in 1:length(tcga_cl_list)){
        #count of samples per tissue type for this cluster
        tissue_types <- tcga_cl_list[[i]]$Tissue_Type %>% unique
        n_per_tissue <- numeric()
        for(j in 1:length(tissue_types)){
            n_per_tissue[j] <- sum(tcga_cl_list[[i]]$Tissue_Type == tissue_types[j])
            names(n_per_tissue)[j] <- tissue_types[j]
        }
        
        # count of samples per group
        # groups <- tcga_cl_list[[i]]$Cluster %>% unique %>% sort
        groups <- clusters_raw[[i]]$Group %>% unique %>% sort
        n_per_group_tissue <- data.frame(
            Tissue_Type = character(), Group = character(), FractionOfTissue = numeric(), Count = numeric())
        this_row <- data.frame(Tissue_Type = 0, Group = 0, FractionOfTissue = 0, Count = 0)
        for(k in 1:length(groups)){
            for(j in 1:length(tissue_types)){
                this_row$Tissue_Type <- tissue_types[j]
                this_row$Group <- groups[k]
                this_row$Count <- tcga_cl_list[[i]][
                    tcga_cl_list[[i]]$Cluster == groups[k] & tcga_cl_list[[i]]$Tissue_Type == tissue_types[j],] %>% nrow
                this_row$FractionOfTissue <- this_row$Count / n_per_tissue[names(n_per_tissue) == tissue_types[j]]
                
                n_per_group_tissue %<>% rbind(this_row) 
            }
            
        }
        
        # plot
        
        tcga_clust_tx_gg_list[[i]] <- ggplot(n_per_group_tissue) +
            geom_bar(aes(x = Tissue_Type, fill = Tissue_Type, y = FractionOfTissue), stat = "identity") +
            facet_grid(.~Group) +
            ggtitle(names(tcga_cl_list)[i]) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(x = Tissue_Type, y = FractionOfTissue, label = Count))
    }
    
    # ----
    # Drug Plot Function
    # ----
    drugPlot <- function(dsdata, pairwise, cluster, clustname, drug){
        #dsdata drug sensitivity table. ds_call_df
        #pairwise : hits or hits_filtered or similar to identify important pairs
        #cluster clusters_raw[[x]] table w/w Cell_Line and Group variables
        #clustname names(clusters_raw)[x] for plot title
        #drug character string drug name in dsdata
        
        pd <- pairwise[pairwise$Drug == drug, ] # subset pairwise data specific to this drug
        sub_ds <- dsdata[dsdata$Drug == drug, ] # subset drug senstivity data specific to this drug
        sub_ds$Group <- cluster$Group[match(sub_ds$Cell_Line, cluster$Cell_Line)] %>% as.factor #annotate ds data with groups from clustering
        
        signif_pairs <- pd$Group_Pair %>% gsub("^G_", "", .) %>% strsplit("-") #process pairwise data to get list of pairs
        
        # For stat_value_manual
        pd$group1 <- signif_pairs %>% lapply(FUN = function(x) x %>% extract(1)) %>% unlist
        pd$group2 <- signif_pairs %>% lapply(FUN = function(x) x %>% extract(2)) %>% unlist
        pd$p <- pd$padj_2 %>% formatC(format = "e", digits = 2)
        pd$y.position <- ceiling(max(sub_ds$LFC)) + (1:nrow(pd))/2
        
        # update any misleading "zero" p values
        nonzero_p <- pd$padj_2[pd$padj_2 > 0]
        min_nonzero_pe <- pd$p[pd$padj_2 == min(nonzero_p)][1]
        pd$p[pd$p == "0.00e+00"] <- paste0("<", min_nonzero_pe)
        
        annotations <- data.frame(
            xpos = c(-Inf,-Inf,Inf,Inf),
            ypos =  c(-Inf, Inf,-Inf,Inf),
        annotateText = c( #bottom left, top left, bottom right, top right
            paste0("Drug MOA: ", dsdata$moa[dsdata$Drug == drug][1]),
            "",
            "", # paste0("PlotID: ", "Clust", clustname, "-", "Tx", drug),
            ""),
        hjustvar = c(0,0,1,1) ,
        vjustvar = c(0,1,0,1)) #<- adjust
        
        txt_size <- 15
        
        g <- ggplot(sub_ds %>% na.omit, aes(x = Group, y = LFC)) +
            
            #Reference Lines
            geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
            geom_hline(yintercept = 1, alpha = 0.5, linetype = "longdash") +
            geom_hline(yintercept = -1, alpha = 0.5, linetype = "longdash") +
            
            
            #Main Data
            geom_boxplot(alpha = 0.5) +
            geom_jitter(width = 0.1, size = 3,
                        aes(color = Tissue_Type)) +
            ggtitle(paste("ClustMethod:", clustname, ", Drug: ", drug)) +
            
            ggpubr::stat_pvalue_manual(data = pd, y.position = "y.position") +
            
            geom_text(data = annotations,
                      aes(x = xpos, y = ypos,
                          hjust = hjustvar, vjust = vjustvar,
                          label = annotateText)) +
            
            ylab("Log2 Fold Change vs. DMSO") +
            xlab("Group ID") +
            
            theme(axis.text = element_text(size = txt_size),
                  axis.title = element_text(size = txt_size),
                  legend.text = element_text(size = txt_size),
                  legend.title = element_text(size = txt_size)
            )
        
        
        
        return(g)
        
        
    }
    
    # ----
    
    # ----
    # Load ANOVA hits data
    # ----
    ds_call_df <- data.table::fread("CCLE/ccle_ds_lfc-plus-binary.csv", data.table = FALSE)
    pairwise_drug_hits <- data.table::fread("CCLE/pairwise_comp_drugs_hits.csv")
    drug_annot <- data.table::fread("CCLE/primary_replicate_collapsed_treatment_info.csv")
    
    simple_drugnames <- ds_call_df$Drug %>% strsplit("_") %>% lapply(function(x) x[1]) %>% unlist
    
    ds_call_df$simple_drugnames <- simple_drugnames
    ds_call_df$moa <- drug_annot$moa[match(simple_drugnames, drug_annot$name)]
    

    
    
    # ----
    # ANOVA based outputs
    # ----
    
    # Choose Minimum Number of Significant Pairwise Comparisons
    hits_per_drug <- reactive({
        hpd_df <- data.frame(drug = pairwise_drug_hits$Drug %>% unique, nhits = 0)
        for(i in 1:nrow(hpd_df)){
            hpd_df$nhits[i] <- nrow(pairwise_drug_hits[pairwise_drug_hits$ClustMethod == input$clustid & Drug == hpd_df$drug[i],])
        }
        hpd_df
    })
    
    max_hits <- reactive({
        max(hits_per_drug()$nhits %>% na.omit)
    })
    
    
    output$choose_n_signif <- renderUI({
        # id max number of hits
        
        sliderInput("n_signif",
                    "Choose Minimum Significant Pairwise Comparisons",
                    min = 1,
                    max = max_hits(),
                    value = 1, 
                    step = 1)
        

    })
    
    
    hits_filtered <- reactive({
        hpd_df_sub <- hits_per_drug()[hits_per_drug()$nhits >= input$n_signif,]
        pairwise_drug_hits[pairwise_drug_hits$Drug %in% hpd_df_sub$drug,]
    })
    
    
    # Choose Drug MOA to Plot
    output$choose_moa <- renderUI({
        

        
        hit_drugs <- 
            hits_filtered()$Drug[hits_filtered()$ClustMethod == input$clustid] %>% unique
        
        hit_drugs_simple <- 
            hit_drugs %>% strsplit("_") %>% lapply(FUN = function(x) x[1]) %>% unlist
        
        hit_drugs_moa <- drug_annot$moa[match(hit_drugs_simple, drug_annot$name)] %>% unique %>% na.omit
        
        selectInput("moa", "Choose MOA:", as.list(hit_drugs_moa))
        
    })
     
    # Choose Specific Drug to Plot
    output$choose_drug <- renderUI({
        moa_drugs <- drug_annot$name[drug_annot$moa == input$moa] %>% unique
        all_hit_drugs <- hits_filtered()$Drug %>% unique
        all_hit_drugs_simple <- 
            all_hit_drugs %>% strsplit("_") %>% lapply(FUN = function(x) x[1]) %>% unlist
        hit_moa_drugs <- all_hit_drugs_simple[all_hit_drugs_simple %in% moa_drugs]
        
        selectInput("drug", "Choose Drug:", as.list(hit_moa_drugs)) 
        
    })
    
    # ----
    # Output
    # ---
    
    output$clustername <- renderText(paste0("Model: ",names(clusters_raw)[input$clustid]))
    output$ccle_barplot <- renderPlot(clust_tx_gg_list[[input$clustid]])
    output$tcga_barplot <- renderPlot(tcga_clust_tx_gg_list[[input$clustid]])
    
    
    
    output$drug_plot <- renderPlot({
        hit_drugs <- 
            hits_filtered()$Drug[hits_filtered()$ClustMethod == input$clustid] %>% unique
        
        hit_drugs_simple <- 
            hit_drugs %>% strsplit("_") %>% lapply(FUN = function(x) x[1]) %>% unlist
        
        complex_drug <- hit_drugs[hit_drugs_simple == input$drug][1]
        
        drugPlot(
            dsdata = ds_call_df,
            pairwise = hits_filtered()[hits_filtered()$ClustMethod == input$clustid],
            cluster = clusters_raw[[input$clustid]],
            clustname = names(clusters_raw)[input$clustid],
            drug = complex_drug
        )
    })
    
})
