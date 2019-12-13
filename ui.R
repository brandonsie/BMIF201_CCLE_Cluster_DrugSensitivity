library(shiny)
library(shinythemes)
library(markdown)

# UI inputs / ouputs
shinyUI(fluidPage(theme = shinytheme("flatly"), 

    # Application title
    titlePanel("CCLE Cluster Drug Sensitivity"),

    # Sidebar 
    tabsetPanel(               
        tabPanel("App", 
                 
            sidebarLayout(
                
                sidebarPanel(
                    
                    # Slider input for clustering model
                    sliderInput("clustid",
                                "Choose Clustering Model",
                                min = 1,
                                max = 7,
                                value = 1),
                    
                    # Model information
                    textOutput("clustername"),
                    textOutput("this_clust_desc"),
                    br(),
                    
                    # Slider input for number of significant pairwise comparisons
                    uiOutput("choose_n_signif"),
                    
                    # Drop down menu for Drug MOA
                    uiOutput("choose_moa"),
                    
                    # Drop down menu for specific drug
                    uiOutput("choose_drug"), 
                    br(),
                    
                    # Info
                    
                    p("A bit more information is available at the ",
                      a("Github Repo", href="https://github.com/brandonsie/BMIF201_CCLE_Cluster_DrugSensitivity")
                    ), 
                    br(),
                    
                    includeMarkdown("markdown/AppInstructions.md"),
                    br()
                    
                ),
                
                mainPanel(
                    h3("CCLE Cluster Drug Sensitivity Plot"),
                    plotOutput("drug_plot"),
                    
                    h3("CCLE Cell Line Cluster Assignments"),
                    plotOutput("ccle_barplot", height = "200px"),
                    h3("TCGA Tumor Sample Cluster Assignments"),
                    plotOutput("tcga_barplot", height = "200px")
                    
                )
            ) # end sidebaryLayout
        ), #end tabPanel
        tabPanel("Information",
                 includeMarkdown("README.md")
        )
    )
))


