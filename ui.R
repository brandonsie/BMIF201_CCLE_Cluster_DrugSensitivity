library(shiny)

# UI inputs / ouputs
shinyUI(fluidPage(

    # Application title
    titlePanel("CCLE Cluster Drug Sensitivity"),

    # Sidebar 
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
    )
))


# 4686 drug / concentration combos, 209 cell lines ( from our 4 tissue types) 
# 367 cell lines in clusters. not all have drug sensitivity. 
# also some cell lines not screened against some drugs. e.g. 207/209 for first drug
