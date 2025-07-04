#!/usr/bin/env Rscript
library(shiny)
library(plotly)
library(data.table)

#command-line args
args <- commandArgs(trailingOnly = TRUE)

#input folder - this is just the cosigt output folder
input_folder <- args[1]
#other inputs are generated by cosigt workflow in specific folders
gafpack_folder <- file.path(input_folder, "gafpack")
odgi_paths_folder <- file.path(input_folder, "odgi/paths")
odgi_len_folder<-file.path(input_folder, "odgi/view")

#list available sample vectors
list_sample_region_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.gafpack\\.gz$", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    sample = basename(dirname(dirname(dirname(files)))),
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  return(file_info)
}

#list available node length tables
list_node_length_files <- function(folder) {
  files <- list.files(folder, pattern = ".node.length.tsv", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  return(file_info)
}

#list available pangenome coverage tables
list_region_tables <- function(folder) {
  files <- list.files(folder, pattern = "\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  return(file_info)
}

#list available mask and submasks files
list_region_mask <- function(folder) {
  files <- list.files(folder, pattern = "mask", full.names = TRUE, recursive = TRUE)
  file_info <- data.frame(
    file_path = files,
    region = do.call(c, lapply(files, function(x) unlist(strsplit(basename(x), ".", fixed=TRUE))[[1]]))
  )
  #adjust
  for (i in 1:nrow(file_info)) {
    if (startsWith(file_info$region[i], "mask")) {
      file_info$region[i]<-unlist(strsplit(basename(dirname(file_info$file_path[i])), "_submasks"))[1]
    }
  }
  return(file_info)
}

#read a gzipped table from file into a vector
load_vector_from_gz <- function(file_path) {
  df <- fread(file_path)
  return(as.numeric(df[, 2:ncol(df)]))
}

#read a gzipped table from file
load_table_from_gz <- function(file_path) {
  return(fread(file_path))
}

#read node lengths from file
load_node_lengths <- function(file_path) {
  node_lengths <- fread(file_path, header = FALSE, col.names = c("node_id", "length"))
  return(node_lengths)
}

#read nodes mask from file
load_mask <- function(file_path) {
  mask <- fread(file_path, header = FALSE, col.names = c("mask"))
  return(mask)
}

#calculate cosine similarity between two vectors
cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

#calculate cumulative lengths for nodes
calculate_cumulative_lengths <- function(lengths) {
  lengths$cumulative_start <- cumsum(c(1, lengths$length[-nrow(lengths)]))
  lengths$cumulative_end <- lengths$cumulative_start + lengths$length
  return(lengths)
}

#load informations from file
sample_file_info <- list_sample_region_files(gafpack_folder)
table_file_info <- list_region_tables(odgi_paths_folder)
node_length_file_info<-list_node_length_files(odgi_len_folder)
mask_file_info<-list_region_mask(odgi_paths_folder)

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

      # Mask file selection for the chosen region
      uiOutput("mask_ui"),
      # Interactive button to toggle mask
      actionButton("toggle_mask", "Toggle Mask"),

      # Row selection for the first plot
      uiOutput("row_select1"),
      uiOutput("row_select2"),

      # Row selection for the second plot
      uiOutput("row_select3"),
      uiOutput("row_select4"),

      # Load and compare button
      actionButton("load_compare", "Load and Compare Vectors"),
      textOutput("cosineSimilarity"),
      textOutput("cosineSimilarity2"),
    ),

    mainPanel(
      plotlyOutput("vectorPlot"),  # First plot
      plotlyOutput("vectorPlot2")  # Second plot
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

  # Filter available mask files based on selected region
  output$mask_ui <- renderUI({
    req(input$region)
    masks <- mask_file_info[mask_file_info$region == input$region, "file_path"]
    mask_names <- basename(masks)
    selectInput("mask", "Mask File:", choices = setNames(masks, mask_names))
  })
  
  # Dynamic UI to select two rows in the chosen table
  output$row_select1 <- renderUI({
    req(input$table)
    table_data <- load_table_from_gz(input$table)
    selectInput("row1", "First Haplotype (upper plot):", choices = table_data$path.name)
  })

  output$row_select2 <- renderUI({
    req(input$table)
    table_data <- load_table_from_gz(input$table)
    selectInput("row2", "Second Haplotype (upper plot):", choices = table_data$path.name)
  })

 # Existing reactive logic for the first plot
  output$row_select3 <- renderUI({
    req(input$table)
    table_data <- load_table_from_gz(input$table)
    selectInput("row3", "First Haplotype (lower plot):", choices = table_data$path.name)
  })

  output$row_select4 <- renderUI({
    req(input$table)
    table_data <- load_table_from_gz(input$table)
    selectInput("row4", "Second Haplotype (lower plot):", choices = table_data$path.name)
  })

  # Reactive values to hold the loaded vectors
  sample_vector <- reactiveVal()
  summed_vector <- reactiveVal()
  pang_vector1 <- reactiveVal()
  pang_vector2 <- reactiveVal()
  pang_vector3 <- reactiveVal()
  pang_vector4 <- reactiveVal()
  summed_vector2 <- reactiveVal()
  node_lengths <- reactiveVal()
  mask_file<-reactiveVal()
  mask_active <- reactiveVal(FALSE)
  
  # Load the mask when the toggle button is clicked
  observeEvent(input$toggle_mask, {
    req(input$mask)
    mask_file(load_mask(input$mask))
    mask_active(!mask_active()) # Toggle the mask AFTER loading the mask
  })

  # Load and compare vectors when button is clicked
  observeEvent(input$load_compare, {
    req(input$sample, input$table, input$row1, input$row2, input$length, input$row3, input$row4)

    # Load the selected sample vector
    sample_file <- sample_file_info$file_path[sample_file_info$sample == input$sample & sample_file_info$region == input$region]
    sample_vector(load_vector_from_gz(sample_file))

    # Load the table and extract the selected rows
    table_data <- load_table_from_gz(input$table)
    row1_vector <- as.numeric(table_data[table_data$path.name == input$row1, -1])
    row2_vector <- as.numeric(table_data[table_data$path.name == input$row2, -1])
    row3_vector <- as.numeric(table_data[table_data$path.name == input$row3, -1])
    row4_vector <- as.numeric(table_data[table_data$path.name == input$row4, -1])

    # Sum the selected rows
    summed_vector(row1_vector + row2_vector)
    summed_vector2(row3_vector + row4_vector)
    #load mask
    mask_file(load_mask(input$mask))

    # Store individual vectors
    pang_vector1(row1_vector)
    pang_vector2(row2_vector)
    pang_vector3(row3_vector)
    pang_vector4(row4_vector)
    lengths <- load_node_lengths(input$length)
    # Calculate cumulative lengths
    lengths <- calculate_cumulative_lengths(lengths)
    node_lengths(lengths)

    # If mask is not loaded yet, load it if available
    if (is.null(mask_file())) {
      req(input$mask)
      mask_file(load_mask(input$mask))
    }

  })

  # Calculate cosine similarity and angle between the sample vector and summed vector
  output$cosineSimilarity <- renderText({
    req(sample_vector(), summed_vector())
    mask <- if (!is.null(mask_file())) mask_file() else data.frame(mask = rep(1, length(sample_vector())))
    apply_mask <- mask_active()

    # Apply the mask
    valid_nodes <- if (apply_mask) mask$mask == 1 else TRUE

    similarity <- cosine_similarity(
      sample_vector()[valid_nodes],
      summed_vector()[valid_nodes]
    )

    angle <- acos(similarity) * (180 / pi)  # Convert to degrees

    HTML(paste("Cosine Similarity (upper plot): ", round(similarity, 4),"; ",
              "Angle (degrees, upper plot): ", round(angle, 2),
              sep = ""))
  })

  output$cosineSimilarity2 <- renderText({
    req(sample_vector(), summed_vector2())
    mask <- if (!is.null(mask_file())) mask_file() else data.frame(mask = rep(1, length(sample_vector())))
    apply_mask <- mask_active()

    # Apply the mask
    valid_nodes <- if (apply_mask) mask$mask == 1 else TRUE

    similarity <- cosine_similarity(
      sample_vector()[valid_nodes],
      summed_vector2()[valid_nodes]
    )

    angle <- acos(similarity) * (180 / pi)  # Convert to degrees

    HTML(paste("Cosine Similarity (lower plot): ", round(similarity, 4),"; ",
              "Angle (degrees, lower plot): ", round(angle, 2),
              sep = ""))
  })

  output$vectorPlot <- renderPlotly({
    req(sample_vector(), summed_vector(), pang_vector1(), pang_vector2(), node_lengths(), mask_file())

    mask <- if (!is.null(mask_file())) mask_file() else data.frame(mask = rep(1, length(sample_vector())))
    apply_mask <- mask_active()

    # Extract node lengths and validate dimensions
    lengths <- node_lengths() 
    node_ids <- seq_along(sample_vector())

    # Mask
    mask <- if (!is.null(mask_file())) mask_file() else data.frame(mask = rep(1, length(sample_vector())))
    apply_mask <- mask_active()
    masked_nodes <- if (apply_mask) mask$mask == 1 else TRUE

    #Normalize by magnitude
    y1_max <- max(sample_vector(), na.rm = TRUE)
    y2_max<- max(summed_vector(), na.rm = TRUE)

    normalize_y1 <- function(x) x / y1_max
    normalize_y2 <- function(x) x / y2_max

    # Prepare the plot data using node midpoints for x-coordinates
    plot_data <- data.frame(
      x = (lengths$cumulative_start + lengths$cumulative_end) / 2,  # Node midpoint
      y = c(sample_vector(),
            summed_vector(),
            pang_vector1(),
            pang_vector2()),
      y_normalized = c(normalize_y1(sample_vector()),
                      normalize_y2(summed_vector()),
                      normalize_y2(pang_vector1()),
                      normalize_y2(pang_vector2())),
      node_id = rep(node_ids, 4),
      Vector = rep(c("Sample Vector (diploid)",
                    "Pangenome Vector (diploid)",
                    "Pangenome Vector (1st hap)",
                    "Pangenome Vector (2nd hap)"),
                  each = length(node_ids)),
      masked = rep(masked_nodes, 4)
    )

    # Create the plot
    plot <- plot_ly()

    # Add trace for Sample Vector (diploid) using primary y-axis (y1)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Sample Vector (diploid)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Sample Vector (diploid)",
        marker = list(size = 5, opacity = ifelse(plot_data$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Sample Vector (diploid)"),
        hoverinfo = "text",
        yaxis = "y1"
      )

    # Add trace for Pangenome Vector (diploid) using secondary y-axis (y2)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Pangenome Vector (diploid)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (diploid)",
        marker = list(size = 5, opacity = ifelse(plot_data$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (diploid)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add trace for Pangenome Vector (1st hap) using secondary y-axis (y2)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Pangenome Vector (1st hap)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (1st hap)",
        marker = list(size = 5, opacity = ifelse(plot_data$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (1st hap)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add trace for Pangenome Vector (2nd hap) using secondary y-axis (y2)
    plot <- plot %>%
      add_trace(
        data = subset(plot_data, Vector == "Pangenome Vector (2nd hap)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (2nd hap)",
        marker = list(size = 5, opacity = ifelse(plot_data$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (2nd hap)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add layout with dual y-axis
    plot <- plot %>%
      layout(
        title = "Comparison of Selected Vectors",
        xaxis = list(title = "Cumulative Haplotype Length"),
        yaxis = list(title = "Sample Vector Coverage (Scaled)", side = "left", range=c(-0.1,1.1)),
        yaxis2 = list(
          title = "Pangenome Vector Coverage (Scaled)",
          side = "right",
          overlaying = "y",  # Overlay on the same x-axis
          showgrid = FALSE,
          range=c(-0.1,1.1)
        ),
        legend = list(x = 1.1, y = 0.5)
      )

    return(plot)
  })

  # Render the second plot
  output$vectorPlot2 <- renderPlotly({
    req(sample_vector(), summed_vector2(), pang_vector3(), pang_vector4(), node_lengths(), mask_file())

    # Extract node lengths and validate dimensions
    lengths <- node_lengths()
    node_ids <- seq_along(sample_vector())

    # Mask
    mask <- if (!is.null(mask_file())) mask_file() else data.frame(mask = rep(1, length(sample_vector())))
    apply_mask <- mask_active()
    masked_nodes <- if (apply_mask) mask$mask == 1 else TRUE

    #Normalize by magnitude
    y1_max <- max(sample_vector(), na.rm = TRUE)
    y2_max<- max(summed_vector2(), na.rm = TRUE)

    normalize_y1 <- function(x) x / y1_max
    normalize_y2 <- function(x) x / y2_max

    # Prepare the plot data for the second plot
    plot_data2 <- data.frame(
      x = (lengths$cumulative_start + lengths$cumulative_end) / 2,  # Node midpoint
      y = c(sample_vector(),
            summed_vector2(),
            pang_vector3(),
            pang_vector4()),
      y_normalized = c(normalize_y1(sample_vector()),
                      normalize_y2(summed_vector2()),
                      normalize_y2(pang_vector3()),
                      normalize_y2(pang_vector4())),
      node_id = rep(node_ids, 4),
      Vector = rep(c("Sample Vector (diploid)",
                    "Pangenome Vector (diploid)",
                    "Pangenome Vector (1st hap)",
                    "Pangenome Vector (2nd hap)"),
                  each = length(node_ids)),
      masked = rep(masked_nodes, 4)
    )

    # Create the plot
    plot2 <- plot_ly()

    # Add trace for Sample Vector (diploid) using primary y-axis (y1)
    plot2 <- plot2 %>%
      add_trace(
        data = subset(plot_data2, Vector == "Sample Vector (diploid)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Sample Vector (diploid)",
        marker = list(size = 5, opacity = ifelse(plot_data2$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Sample Vector (diploid)"),
        hoverinfo = "text",
        yaxis = "y1"
      )

    # Add trace for Pangenome Vector (diploid, Plot 2)
    plot2 <- plot2 %>%
      add_trace(
        data = subset(plot_data2, Vector == "Pangenome Vector (diploid)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (diploid)",
        marker = list(size = 5, opacity = ifelse(plot_data2$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (diploid)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add trace for Pangenome Vector (1st hap, Plot 2)
    plot2 <- plot2 %>%
      add_trace(
        data = subset(plot_data2, Vector == "Pangenome Vector (1st hap)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (1st hap)",
        marker = list(size = 5, opacity = ifelse(plot_data2$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (1st hap)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add trace for Pangenome Vector (2nd hap, Plot 2)
    plot2 <- plot2 %>%
      add_trace(
        data = subset(plot_data2, Vector == "Pangenome Vector (2nd hap)"),
        x = ~x, y = ~y_normalized,
        type = "scatter",
        mode = "markers",
        name = "Pangenome Vector (2nd hap)",
        marker = list(size = 5, opacity = ifelse(plot_data2$masked, 0.9, 0.2)),
        text = ~paste("Node ID:", node_id, "<br>Coverage:", y, "<br>Vector: Pangenome Vector (2nd hap)"),
        hoverinfo = "text",
        yaxis = "y2"
      )

    # Add layout with dual y-axis
    plot2 <- plot2 %>%
      layout(
        xaxis = list(title = "Cumulative Haplotype Length"),
        yaxis = list(title = "Sample Vector Coverage (Scaled)", side = "left", range=c(-0.1,1.1)),
        yaxis2 = list(
          title = "Pangenome Vector Coverage (Scaled)",
          side = "right",
          overlaying = "y",  # Overlay on the same x-axis
          showgrid = FALSE,
          range=c(-0.1,1.1)
        ),
        legend = list(x = 1.1, y = 0.5)
      )

    return(plot2)
  })

}

# Run the Shiny app
shinyApp(ui = ui, server = server)