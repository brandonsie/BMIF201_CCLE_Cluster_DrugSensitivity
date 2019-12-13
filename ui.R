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
            br(),
            
            uiOutput("choose_n_signif"),
            
            uiOutput("choose_moa"),
            uiOutput("choose_drug"),
            
            br(),
            
            
            br()
            
        ),

        
        
        # Show a plot of the generated distribution
        mainPanel(
            
            #plotOutput("distPlot")
            # tableOutput("dimThis")
            plotOutput("drug_plot"),
            
            p("CCLE Cell Line Cluster Assignments"),
            plotOutput("ccle_barplot"),
            p("TCGA Tumor Sample Cluster Assignments"),
            plotOutput("tcga_barplot"),
            
            tableOutput("hits_filtered_head"),
            tableOutput("nsignif")
            
        )
    )
))


# 4686 drug / concentration combos, 209 cell lines ( from our 4 tissue types) 
# 367 cell lines in clusters. not all have drug sensitivity. 
# also some cell lines not screened against some drugs. e.g. 207/209 for first drug
