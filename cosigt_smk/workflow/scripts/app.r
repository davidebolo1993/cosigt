library(shiny)
library(plotly)
library(data.table)
library(stringr)  # For handling strings

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input folder containing gafpack vectors
input_folder <- args[1]
odgi_paths_folder <- file.path(input_folder, "odgi/paths/matrix")

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

# Function to calculate cosine similarity
cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

# Load file information for vectors and tables
sample_file_info <- list_sample_region_files(input_folder)
table_file_info <- list_region_tables(odgi_paths_folder)

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

  # Load and compare vectors when button is clicked
  observeEvent(input$load_compare, {
    req(input$sample, input$table, input$row1, input$row2)

    # Load the selected sample vector
    sample_file <- sample_file_info$file_path[sample_file_info$sample == input$sample & sample_file_info$region == input$region]
    sample_vector(load_vector_from_gz(sample_file))

    # Load the table and extract the selected rows
    table_data <- load_table_from_gz(input$table)
    row1_vector <- as.numeric(table_data[table_data$path.name == input$row1, -1])
    row2_vector <- as.numeric(table_data[table_data$path.name == input$row2, -1])

    # Sum the selected rows
    summed_vector(row1_vector + row2_vector)
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

  # Plot the two selected vectors
  output$vectorPlot <- renderPlotly({
    req(sample_vector(), summed_vector())

    # Data for plotting
    plot_data <- data.frame(
      x = c(1:length(sample_vector()), 1:length(summed_vector())),
      y = c(sample_vector(), summed_vector()),
      Vector = rep(c("Sample Vector", "Pangenome Vector"), each = length(sample_vector()))
    )

    # Plot using plotly
    plot <- plot_ly(plot_data, x = ~x, y = ~y, color = ~Vector, type = "scatter", mode = "lines+markers") %>%
      layout(title = paste("Comparison of Selected Vectors"),
             xaxis = list(title = "Index"),
             yaxis = list(title = "Value"))

    plot
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

