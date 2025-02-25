# Load necessary libraries
library(AdhereR)
library(dplyr)
library(plyr)
library(lubridate)
library(latticeExtra)
library(data.table)
library(factoextra)
library(stats)

# Prepare the dataset
ExamplePats <- med.events
tidy_data <- ExamplePats
colnames(tidy_data) <- c("patient_id", "event_date", "dose_per_day", "ATC_code", "original_duration")
tidy_data$event_date <- mdy(tidy_data$event_date)

# Function to analyze medication events
analyze_medication <- function(medication_code) {
  # Filter data for the specified medication
  medication_data <- tidy_data %>% filter(ATC_code == medication_code)
  
  # Arrange and group data by patient and event date
  processed_data <- medication_data %>%
    arrange(patient_id, event_date) %>%
    group_by(patient_id) %>%
    mutate(previous_event_date = lag(event_date, default = NA))
  
  # Remove rows with NA in previous_event_date
  processed_data <- processed_data %>% filter(!is.na(previous_event_date))
  
  # Sample one random row per patient
  sampled_data <- ddply(processed_data, .(patient_id), function(x) x[sample(nrow(x), 1),])
  sampled_data <- sampled_data %>% select(patient_id, event_date, previous_event_date)
  
  # Calculate event intervals
  sampled_data$event_interval <- as.numeric(sampled_data$event_date - sampled_data$previous_event_date)
  
  # Generate ECDF plot
  ecdf_plot <- ecdfplot(~sampled_data$event_interval)
  ecdf_values <- ecdf_plot$panel.args[[1]]
  ecdf_functions <- lapply(split(sampled_data$event_interval, 1), ecdf)
  ecdf_results <- sapply(ecdf_functions, function(e) e(sampled_data$event_interval))
  ecdf_results <- as.vector(ecdf_results)
  ecdf_values <- unlist(ecdf_values)
  ecdf_values <- as.numeric(ecdf_values)
  
  # Create a data frame for ECDF results
  ecdf_df <- data.frame(x = ecdf_values, y = ecdf_results)
  
  # Retain the lower 80% of the ECDF
  ecdf_df <- ecdf_df %>% filter(y <= 0.8)
  
  # Plot ECDFs
  par(mfrow = c(1, 2))
  plot(ecdf_df$x, ecdf_df$y, main = "80% ECDF")
  plot(ecdf_values, ecdf_results, main = "100% ECDF")
  
  # Calculate maximum interval in the 80% ECDF
  max_interval <- max(ecdf_df$x)
  
  # Filter data based on the maximum interval
  filtered_data <- sampled_data %>% filter(event_interval <= max_interval)
  
  # Density plot of log-transformed event intervals
  density_plot <- density(log(filtered_data$event_interval))
  plot(density_plot, main = "Log(event interval)")
  
  # Prepare data for clustering
  scaled_data <- scale(data.table(x = density_plot$x, y = density_plot$y))
  
  # Silhouette analysis to determine optimal clusters
  set.seed(1234)
  silhouette_plot <- fviz_nbclust(scaled_data, kmeans, method = "silhouette") + labs(subtitle = "Silhouette Analysis")
  plot(silhouette_plot)
  
  # Determine the optimal number of clusters
  optimal_clusters <- as.numeric(silhouette_plot$data$clusters[which.max(silhouette_plot$data$y)])
  
  # Perform K-means clustering
  set.seed(1234)
  kmeans_result <- kmeans(ecdf_df$x, optimal_clusters)
  ecdf_df$cluster <- as.numeric(kmeans_result$cluster)
  
  # Summarize cluster results
  cluster_summary <- data.frame(
    Cluster = names(tapply(log(ecdf_df$x), ecdf_df$cluster, min)),
    Minimum = exp(tapply(log(ecdf_df$x), ecdf_df$cluster, min)),
    Maximum = exp(tapply(log(ecdf_df$x), ecdf_df$cluster, max)),
    Median = exp(tapply(log(ecdf_df$x), ecdf_df$cluster, median, na.rm = TRUE))
  )
  
  # Merge cluster results with the original data
  final_results <- sampled_data %>%
    cross_join(cluster_summary) %>%
    mutate(Final_cluster = ifelse(event_interval >= Minimum & event_interval <= Maximum, Cluster, NA))
  
  final_results <- final_results %>% filter(!is.na(Final_cluster))
  final_results <- final_results %>% select(patient_id, Median, Cluster)
  
  # Determine the most frequent cluster
  most_frequent_cluster <- as.data.frame(table(final_results$Cluster)) %>%
    arrange(desc(Freq)) %>%
    slice(1) %>%
    select(Var1) %>%
    rename(Cluster = Var1)
  
  # Merge most frequent cluster with final results
  final_results <- merge(most_frequent_cluster, final_results, by = "Cluster")
  final_results <- final_results %>% slice(1) %>% select(-Cluster)
  
  # Assign median and cluster to the original data
  processed_data <- merge(processed_data, final_results, by = "patient_id", all.x = TRUE)
  processed_data$Median <- ifelse(is.na(processed_data$Median), final_results$Median, processed_data$Median)
  processed_data$Cluster <- ifelse(is.na(processed_data$Cluster), "0", processed_data$Cluster)
  processed_data$event_interval <- as.numeric(processed_data$event_interval)
  processed_data$test <- round(processed_data$event_interval - processed_data$Median, 1)
  
  # Return the processed data
  return(processed_data)
}

# Function to visualize medication assumptions
visualize_assumptions <- function(medication_data) {
  medication_data <- medication_data %>%
    arrange(patient_id, event_date) %>%
    group_by(patient_id) %>%
    mutate(previous_event_date = lag(event_date, default = NA))
  
  medication_data <- medication_data %>%
    group_by(patient_id) %>%
    arrange(patient_id, event_date) %>%
    mutate(prescription_number = seq_along(event_date))
  
  medication_data <- medication_data %>% filter(prescription_number >= 2)
  medication_data <- medication_data %>% select(patient_id, event_date, previous_event_date, prescription_number)
  medication_data$Duration <- as.numeric(medication_data$event_date - medication_data$previous_event_date)
  medication_data$prescription_number <- as.factor(medication_data$prescription_number)
  
  # Calculate median duration per patient
  median_durations <- medication_data %>%
    group_by(patient_id) %>%
    summarise(median_duration = median(Duration, na.rm = TRUE))
  
  # Plot boxplot with median durations
  plot <- ggplot(medication_data, aes(x = prescription_number, y = Duration)) +
    geom_boxplot() +
    geom_hline(yintercept = median_durations$median_duration, linetype = "dashed", color = "red") +
    theme_bw()
  
  return(plot)
}

# Analyze medications
medA_data <- analyze_medication("medA")
medB_data <- analyze_medication("medB")

# Visualize assumptions
visualize_assumptions(medA_data)
visualize_assumptions(medB_data)
