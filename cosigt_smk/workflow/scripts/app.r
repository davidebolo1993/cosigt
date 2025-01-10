#!/usr/bin/env Rscript

library(shiny)
library(plotly)
library(data.table)
library(stringr)  # For handling strings

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input folder containing gafpack vectors
input_folder <- args[1]
odgi_paths_folder <- file.path(input_folder, "odgi/paths/matrix")
odgi_len_folder<-file.path(input_folder, "odgi/view")

# Function to list available sample vector files with sample and region info
list_sample_region_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.gafpack\\.gz$", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    sample = basename(dirname(files)),
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  return(file_info)
}

# Function to list available node length tables
list_node_length_files <- function(folder) {
  files <- list.files(folder, pattern = ".node.length.tsv", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  return(file_info)
}


# Function to list available tables with region info
list_region_tables <- function(folder) {
  files <- list.files(folder, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  return(file_info)
}

# Function to read a gzipped vector file
load_vector_from_gz <- function(file_path) {
  df <- fread(file_path)
  return(as.numeric(df[, 2:ncol(df)]))
}

# Function to read a gzipped table
load_table_from_gz <- function(file_path) {
  return(fread(file_path))
}

# Function to read node lengths from file
load_node_lengths <- function(file_path) {
  node_lengths <- fread(file_path, header = FALSE, col.names = c("node_id", "length"))
  return(node_lengths)
}

# Function to calculate cosine similarity
cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

# Function to calculate cumulative lengths for nodes
calculate_cumulative_lengths <- function(lengths) {
  lengths$cumulative_start <- cumsum(c(1, lengths$length[-nrow(lengths)]))
  lengths$cumulative_end <- lengths$cumulative_start + lengths$length
  return(lengths)
}

# Load file information for vectors and tables
sample_file_info <- list_sample_region_files(input_folder)
table_file_info <- list_region_tables(odgi_paths_folder)
node_length_file_info<-list_node_length_files(odgi_len_folder)

# Shiny UI
ui <- fluidPage(
  titlePanel("Read coverage vs Pangenome coverage visualisation"),

  sidebarLayout(
    sidebarPanel(
      # Region selection (applies to both sample vectors and tables)
      selectInput("region", "Region:", choices = unique(sample_file_info$region)),

      # Sample selection for the chosen region
      uiOutput("sample_ui"),

      # Table selection for the chosen region
      uiOutput("table_ui"),

      # Node length file selection for the chosen region
      uiOutput("length_ui"),

      # Row selection within the chosen table (for vector summing)
      uiOutput("row_select1"),
      uiOutput("row_select2"),

      # Load and compare button
      actionButton("load_compare", "Load and Compare Vectors"),
      textOutput("cosineSimilarity"),
      textOutput("angle")
    ),

    mainPanel(
      plotlyOutput("vectorPlot")  # Correctly refer to vectorPlot here
    )
  )
)

# Shiny Server
server <- function(input, output, session) {

  # Filter sample vectors based on selected region
  output$sample_ui <- renderUI({
    req(input$region)
    samples <- sample_file_info[sample_file_info$region == input$region, "sample"]
    selectInput("sample", "Sample:", choices = samples)
  })

  # Filter available tables based on selected region
  output$table_ui <- renderUI({
    req(input$region)
    tables <- table_file_info[table_file_info$region == input$region, "file_path"]
    table_names <- basename(tables)
    selectInput("table", "Pangenome:", choices = setNames(tables, table_names))
  })

  # Filter available node length files based on selected region
  output$length_ui <- renderUI({
    req(input$region)
    lengths <- node_length_file_info[node_length_file_info$region == input$region, "file_path"]
    length_names <- basename(lengths)
    selectInput("length", "Node Length File:", choices = setNames(lengths, length_names))
  })

  # Dynamic UI to select two rows in the chosen table
  output$row_select1 <- renderUI({
    req(input$table)
    table_data <- load_table_from_gz(input$table)
    selectInput("row1", "First Haplotype:", choices = table_data$path.name)
  })

  output$row_select2 <- renderUI({
    req(input$table)
    table_data <- load_table_from_gz(input$table)
    selectInput("row2", "Second Haplotype:", choices = table_data$path.name)
  })

  # Reactive values to hold the loaded vectors
  sample_vector <- reactiveVal()
  summed_vector <- reactiveVal()
  pang_vector1 <- reactiveVal()
  pang_vector2 <- reactiveVal()
  # Reactive value to hold node lengths
  node_lengths <- reactiveVal()

  # Load and compare vectors when button is clicked
  observeEvent(input$load_compare, {
    req(input$sample, input$table, input$row1, input$row2, input$length)

    # Load the selected sample vector
    sample_file <- sample_file_info$file_path[sample_file_info$sample == input$sample & sample_file_info$region == input$region]
    sample_vector(load_vector_from_gz(sample_file))

    # Load the table and extract the selected rows
    table_data <- load_table_from_gz(input$table)
    row1_vector <- as.numeric(table_data[table_data$path.name == input$row1, -1])
    row2_vector <- as.numeric(table_data[table_data$path.name == input$row2, -1])

    # Sum the selected rows
    summed_vector(row1_vector + row2_vector)
    
    # Store individual vectors
    pang_vector1(row1_vector)
    pang_vector2(row2_vector)

    lengths <- load_node_lengths(input$length)
  
    # Calculate cumulative lengths
    lengths <- calculate_cumulative_lengths(lengths)
    node_lengths(lengths)

  })

  # Calculate cosine similarity and angle between the sample vector and summed vector
  output$cosineSimilarity <- renderText({
    req(sample_vector(), summed_vector())
    similarity <- cosine_similarity(sample_vector(), summed_vector())
    paste("Cosine Similarity:", round(similarity, 4))
  })

  output$angle <- renderText({
    req(sample_vector(), summed_vector())
    similarity <- cosine_similarity(sample_vector(), summed_vector())
    angle <- acos(similarity) * (180 / pi)  # Convert to degrees
    paste("Angle (degrees):", round(angle, 2))
  })

  output$vectorPlot <- renderPlotly({
    req(sample_vector(), summed_vector(), pang_vector1(), pang_vector2(), node_lengths())

    # Extract node lengths and validate dimensions
    lengths <- node_lengths()
    node_ids <- seq_along(sample_vector())

    # Ensure lengths and vectors align
    if (nrow(lengths) < length(node_ids)) {
      stop("Node lengths do not match the length of the sample vector.")
    }

    # Prepare the plot data using node midpoints for x-coordinates
    plot_data <- data.frame(
      x = (lengths$cumulative_start + lengths$cumulative_end) / 2,  # Node midpoint
      y = c(sample_vector(),
            summed_vector(),
            pang_vector1(),
            pang_vector2()),
      node_id = rep(node_ids, 4),
      Vector = rep(c("Sample Vector (diploid)",
                    "Pangenome Vector (diploid)",
                    "Pangenome Vector (1st hap)",
                    "Pangenome Vector (2nd hap)"),
                  each = length(node_ids))
    )

    # Create the plot
    plot <- plot_ly()

    # Add trace for Sample Vector (diploid) using primary y-axis (y1)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Sample Vector (diploid)"),
        x = ~x, y = ~y,
        type = "scatter",
        mode = "markers",
        name = "Sample Vector (diploid)",
        marker = list(size = 5, opacity = 0.8),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Sample Vector (diploid)"),
        hoverinfo = "text",
        yaxis = "y1"
      )

    # Add trace for Pangenome Vector (diploid) using secondary y-axis (y2)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Pangenome Vector (diploid)"),
        x = ~x, y = ~y,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (diploid)",
        marker = list(size = 5, opacity = 0.8),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (diploid)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add trace for Pangenome Vector (1st hap) using secondary y-axis (y2)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Pangenome Vector (1st hap)"),
        x = ~x, y = ~y,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (1st hap)",
        marker = list(size = 5, opacity = 0.8),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (1st hap)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add trace for Pangenome Vector (2nd hap) using secondary y-axis (y2)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Pangenome Vector (2nd hap)"),
        x = ~x, y = ~y,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (2nd hap)",
        marker = list(size = 5, opacity = 0.8),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (2nd hap)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add layout with dual y-axis
    plot <- plot %>%
      layout(
        title = "Comparison of Selected Vectors",
        xaxis = list(title = "Cumulative Haplotype Length"),
        yaxis = list(title = "Sample Vector Coverage", side = "left"),
        yaxis2 = list(
          title = "Pangenome Vector Coverage",
          side = "right",
          overlaying = "y",  # Overlay on the same x-axis
          showgrid = FALSE
        ),
        legend = list(x = 1.1, y = 0.5)
      )

    return(plot)
  })


}

# Run the Shiny app
shinyApp(ui = ui, server = server)
