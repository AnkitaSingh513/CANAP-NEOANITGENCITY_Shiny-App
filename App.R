# Load Libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Set Working Directory
setwd("/home/ankitasingh/Desktop/NEOANITGENCITY")

# Load Data
data <- fread("merged_cancer_gene_mutation_info_with_GO.csv")

# Define a function for dynamic color palette
generate_colors <- function(n) {
  if (n <= 8) {
    return(RColorBrewer::brewer.pal(n, "Set2"))
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n))
  }
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = div(
      tags$img(src = "logo.png", height = "50px", style = "margin-right: 10px;"),
      "CANAP - Cancer Neo-Antigen Prediction"
    ),
    titleWidth = 500
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("info-circle")),
      menuItem("Filters & Visualizations", tabName = "filters_visuals", icon = icon("filter")),
      menuItem("Scatter Plot", tabName = "scatter", icon = icon("chart-line")),
      menuItem("GO Analysis", tabName = "go_analysis", icon = icon("chart-bar")),
      menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
      menuItem("Cancer Comparison", tabName = "cancer_comparison", icon = icon("exchange-alt")),
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("Contact Us", tabName = "contact", icon = icon("envelope"))
    )
  ),
  dashboardBody(
    tabItems(
      # Overview Section
      tabItem(tabName = "overview",
              fluidRow(
                box(title = "Overview", width = 12, status = "primary",
                    p("Welcome to CANAP - Cancer Neo-Antigen Prediction Platform!"),
                    p("Analyze cancer gene mutations, visualize mutation frequencies, explore Gene Ontology (GO) terms, 
                       and compare mutation distributions across different cancer types."),
                    p("Key Features:"),
                    tags$ul(
                      tags$li("Filter and visualize gene mutation data with interactive plots."),
                      tags$li("Explore scatter plots for mutation frequency and log difference analysis."),
                      tags$li("Perform GO-based analysis and heatmap visualization."),
                      tags$li("Compare top genes across different cancer types."),
                      tags$li("Download filtered datasets for further analysis.")
                    ),
                    p("This platform supports research in neo-antigenicity and mutation-related analysis.")
                )
              ),
              fluidRow(
                box(title = "Abstract Figure", width = 12, status = "info",
                    tags$img(src = "abstract_figure.jpeg", alt = "Abstract Figure", style = "width:100%;height:auto;")
                )
              )
      ),
      
      # Filters & Visualizations
      tabItem(tabName = "filters_visuals",
              fluidRow(
                box(title = "Filter Options", width = 12, status = "primary",
                    selectizeInput("gene_filter", "Select Gene(s)", choices = NULL, multiple = TRUE),
                    selectizeInput("cancer_filter", "Select Cancer Type(s)", choices = NULL, multiple = TRUE),
                    numericInput("mutation_min", "Min Mutations", value = 0, min = 0, step = 1),
                    numericInput("mutation_max", "Max Mutations", value = 1000, min = 0, step = 1),
                    numericInput("log_diff_min", "Min Log Difference", value = 0, min = 0, step = 0.1),
                    numericInput("log_diff_max", "Max Log Difference", value = 10, min = 0, step = 0.1),
                    selectizeInput("go_filter", "Select GO Term(s)", choices = NULL, multiple = TRUE)
                )
              ),
              fluidRow(
                box(title = "Filtered Data Table", width = 12, status = "primary",
                    DTOutput("filtered_data_table"))
              ),
              fluidRow(
                box(title = "Mutation Frequency (Bar Plot)", width = 6, status = "primary",
                    plotlyOutput("mutation_bar_plot")),
                box(title = "Cancer Type Distribution (Pie Chart)", width = 6, status = "primary",
                    plotlyOutput("cancer_type_plot"))
              )),
      
      # Scatter Plot
      tabItem(tabName = "scatter",
              fluidRow(
                box(title = "Scatter Plot for Gene Analysis", width = 15, status = "primary",
                    plotlyOutput("scatter_plot"))
              )),
      
      # GO Analysis
      tabItem(tabName = "go_analysis",
              fluidRow(
                box(title = "Top 5 GO Categories", width = 6, status = "primary",
                    plotlyOutput("go_bar_plot")),
                box(title = "GO-Based Heatmap", width = 6, status = "primary",
                    plotlyOutput("go_heatmap_plot"))
              )),
      
      # Heatmap
      tabItem(tabName = "heatmap",
              fluidRow(
                box(title = "Top Genes Mutation Heatmap", width = 12, status = "primary",
                    numericInput("top_genes", "Number of Top Genes", value = 10, min = 5, max = 50, step = 5),
                    plotlyOutput("mutation_heatmap_plot"))
              )),
      
      # Cancer Comparison
      tabItem(tabName = "cancer_comparison",
              fluidRow(
                box(title = "Top 5 Genes Across Cancers", width = 12, status = "primary",
                    plotlyOutput("cancer_comparison_plot"))
              )),
      
      # Contact Us Section
      tabItem(tabName = "contact",
              fluidRow(
                box(title = "Contact Us", width = 12, status = "primary",
                    p("For inquiries, contact:"),
                    p(tags$b("Ankita Singh")),
                    p("Email: ", tags$a(href = "mailto:ankita.singh@tii.ae", "ankita.singh@tii.ae")),
                    p(tags$b("Raghvendra Mall")),
                    p("Email: ", tags$a(href = "mailto:raghvendra.mall@tii.ae", "raghvendra.mall@tii.ae"))
                )
              ))
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Update Selectize Inputs Dynamically
  updateSelectizeInput(session, "gene_filter", choices = unique(data$Hugo_Symbol), server = TRUE)
  updateSelectizeInput(session, "cancer_filter", choices = unique(data$project_id), server = TRUE)
  updateSelectizeInput(session, "go_filter", choices = unique(unlist(strsplit(data$GO_Cellular_Component, " & "))), server = TRUE)
  
  # Reactive Data
  filtered_data <- reactive({
    df <- data
    if (!is.null(input$gene_filter)) df <- df[df$Hugo_Symbol %in% input$gene_filter, ]
    if (!is.null(input$cancer_filter)) df <- df[df$project_id %in% input$cancer_filter, ]
    df <- df[df$altered >= input$mutation_min & df$altered <= input$mutation_max, ]
    df <- df[df$log_difference >= input$log_diff_min & df$log_difference <= input$log_diff_max, ]
    if (!is.null(input$go_filter)) {
      df <- df[sapply(strsplit(df$GO_Cellular_Component, " & "), function(x) any(x %in% input$go_filter)), ]
    }
    validate(need(nrow(df) > 0, "No data available for the selected filters."))
    return(df)
  })
  
  # Render Data Table
  output$filtered_data_table <- renderDT({
    datatable(filtered_data(), options = list(scrollX = TRUE))
  })
  
  # Mutation Frequency Plot
  output$mutation_bar_plot <- renderPlotly({
    df <- filtered_data()
    req(nrow(df) > 0)
    plot_ly(df, x = ~Hugo_Symbol, y = ~altered, type = "bar", color = ~project_id,
            colors = generate_colors(length(unique(df$project_id)))) %>%
      layout(title = "Mutation Frequency per Gene", xaxis = list(title = "Genes"), yaxis = list(title = "Mutations"))
  })
  
  # Cancer Type Distribution
  output$cancer_type_plot <- renderPlotly({
    df <- filtered_data()
    req(nrow(df) > 0)
    cancer_summary <- df %>%
      group_by(project_id) %>%
      summarise(Count = n())
    plot_ly(cancer_summary, labels = ~project_id, values = ~Count, type = "pie",
            colors = generate_colors(length(unique(cancer_summary$project_id)))) %>%
      layout(title = "Cancer Type Distribution")
  })
  
  # GO Analysis: Top 5 GO Categories
  output$go_bar_plot <- renderPlotly({
    df <- filtered_data()
    req(nrow(df) > 0)
    go_summary <- df %>%
      separate_rows(GO_Cellular_Component, sep = " & ") %>%
      count(GO_Cellular_Component, name = "Count") %>%
      arrange(desc(Count)) %>%
      slice_head(n = 5)
    plot_ly(go_summary, x = ~GO_Cellular_Component, y = ~Count, type = "bar", color = ~GO_Cellular_Component) %>%
      layout(title = "Top 5 GO Categories", xaxis = list(title = "GO Categories"), yaxis = list(title = "Gene Count"))
  })
  
  # GO Analysis: GO-Based Heatmap
  output$go_heatmap_plot <- renderPlotly({
    df <- filtered_data()
    req(nrow(df) > 0)  # Ensure filtered data exists
    
    # Extract and reshape data for the heatmap
    go_heatmap <- df %>%
      separate_rows(GO_Cellular_Component, sep = " & ") %>%  # Split GO terms
      count(Hugo_Symbol, GO_Cellular_Component) %>%         # Count occurrences
      pivot_wider(names_from = GO_Cellular_Component,       # Reshape to wide format
                  values_from = n, 
                  values_fill = 0)
    
    # Convert to matrix for Plotly heatmap
    heatmap_matrix <- as.matrix(go_heatmap[, -1])  # Exclude the gene column
    rownames(heatmap_matrix) <- go_heatmap$Hugo_Symbol  # Set row names
    
    # Ensure the matrix has data
    validate(need(ncol(heatmap_matrix) > 0 && nrow(heatmap_matrix) > 0, 
                  "No data available for the GO heatmap."))
    
    # Plot the heatmap
    plot_ly(
      z = heatmap_matrix, 
      x = colnames(heatmap_matrix), 
      y = rownames(heatmap_matrix), 
      type = "heatmap", 
      colorscale = "Viridis"
    ) %>%
      layout(
        title = "GO-Based Heatmap",
        xaxis = list(title = "GO Categories"),
        yaxis = list(title = "Genes")
      )
  })
  
  
    # Cancer Comparison Plot
    output$cancer_comparison_plot <- renderPlotly({
      df <- filtered_data()
      req(nrow(df) > 0)
      top_genes <- df %>%
        group_by(project_id, Hugo_Symbol) %>%
        summarise(Total_Mutations = sum(altered, na.rm = TRUE)) %>%
        arrange(desc(Total_Mutations)) %>%
        group_by(project_id) %>%
        slice_head(n = 5)
      plot_ly(top_genes, 
              x = ~Hugo_Symbol, 
              y = ~Total_Mutations, 
              type = "bar", 
              color = ~project_id, 
              colors = generate_colors(length(unique(top_genes$project_id)))) %>%
        layout(title = "Top 5 Genes Across Cancers",
               xaxis = list(title = "Genes"),
               yaxis = list(title = "Total Mutations"),
               barmode = "group")
    })
    
    # Scatter Plot
    output$scatter_plot <- renderPlotly({
      df <- filtered_data()
      req(nrow(df) > 0)
      plot_ly(df, 
              x = ~log_difference, 
              y = ~altered, 
              type = "scatter", 
              mode = "markers",
              text = ~paste("Gene:", Hugo_Symbol, "<br>Mutations:", altered),
              color = ~Hugo_Symbol,
              colors = generate_colors(length(unique(df$Hugo_Symbol)))) %>%
        layout(title = "Scatter Plot: Log Difference vs Mutations",
               xaxis = list(title = "Log Difference"), 
               yaxis = list(title = "Number of Mutations"))
    })
    
    # Heatmap: Mutation Heatmap for Top Genes
    output$mutation_heatmap_plot <- renderPlotly({
      df <- filtered_data()
      req(nrow(df) > 0)
      top_genes <- df %>%
        arrange(desc(altered)) %>%
        head(input$top_genes)
      mutation_heatmap <- top_genes %>%
        select(Hugo_Symbol, Missense_Mutation, Frame_Shift_Del, Nonsense_Mutation, Splice_Site) %>%
        pivot_longer(cols = -Hugo_Symbol, names_to = "Mutation_Type", values_to = "Count") %>%
        acast(Hugo_Symbol ~ Mutation_Type, value.var = "Count", fill = 0)
      plot_ly(z = as.matrix(mutation_heatmap), 
              x = colnames(mutation_heatmap), 
              y = rownames(mutation_heatmap), 
              type = "heatmap", 
              colorscale = "Blues") %>%
        layout(title = "Top Genes Mutation Heatmap", 
               xaxis = list(title = "Mutation Types"), 
               yaxis = list(title = "Genes"))
    })
    
    # Download Data
    output$download_filtered <- downloadHandler(
      filename = "filtered_data.csv",
      content = function(file) {
        write.csv(filtered_data(), file, row.names = FALSE)
      }
    )
}

# Run App
shinyApp(ui = ui, server = server)

            