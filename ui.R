library(shiny)
library(shinythemes)

# UI inputs / ouputs
shinyUI(fluidPage(theme = shinytheme("flatly"), 

    # Application title
    titlePanel("CCLE Cluster Drug Sensitivity"),

    # Sidebar 
    tabsetPanel(               
        tabPanel("App", 
                 
            sidebarLayout(
                
                # lider input for number of bins
                sidebarPanel(
                    sliderInput("clustid",
                                "Choose Clustering Model",
                                min = 1,
                                max = 7,
                                value = 1),
                    
                    textOutput("clustername"),
                    textOutput("this_clust_desc"),
                    br(),
                    
                    uiOutput("choose_n_signif"),
                    
                    uiOutput("choose_moa"),
                    uiOutput("choose_drug"),
                    
                    br(),
                    
                    # Info
                    
                    p(
                        "A bit more information is available at the ",
                        a("Github Repo", href="https://github.com/brandonsie/BMIF201_CCLE_Cluster_DrugSensitivity")
                    ),
                    br()
                    
                ),
        
                
                
                # Show a plot of the generated distribution
                mainPanel(
                    
                    #plotOutput("distPlot")
                    # tableOutput("dimThis")
                    h3("CCLE Cluster Drug Sensitivity Plot"),
                    plotOutput("drug_plot"),
                    
                    h3("CCLE Cell Line Cluster Assignments"),
                    plotOutput("ccle_barplot", height = "200px"),
                    h3("TCGA Tumor Sample Cluster Assignments"),
                    plotOutput("tcga_barplot", height = "200px")
                    
                    
                )
            ) # end sidebaryLayout
        ), #end tabPanel
        tabPanel("Documentation",
                 p(
                     "Extensive credit to Maha Shady, Greg Brunette, Katherine Duchinski", br(), 
                     "One-way ANOVA run on PRISM drug screen replicate collapsed log2 fold change vs. DMSO data for each cluster model",
                     "(4686 drugs, 209 CCLE cell lines (from our 4 tissue types).", br(),
                     "Tukey pairwise comparisons and BH multiple hypothesis correction performed.", br(),
                     "Each ggplot figure represents the drug sensitivity for one drug across the clusters of one model. A figure can be generated if at least one pair of clusters has an adjusted p value < 0.05."
                     
                 ),
                 
                 p(
                     "A bit more information is available at the ",
                     a("Github Repo", href="https://github.com/brandonsie/BMIF201_CCLE_Cluster_DrugSensitivity")
                 )
        )
    )
))


