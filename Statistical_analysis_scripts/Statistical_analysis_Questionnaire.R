# Required libraries
library(readxl)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

library(tidyverse)
library(rstatix)

# PLOTTING QUESTIONNAIRE ENVIRONMENTS ANSWERS ----

# Function to process and plot data from an .xlsx file
questionnaire_environment_plot_xlsx <- function(file_path, sheet, plot_title) {
  
  #data <- read_excel(file_path)
  data <- read_excel(file_path, sheet = sheet)
  
  h_msk <- data$'Highway #1'
  h_spb <- data$'Highway #2'
  b_msk <- data$'Boulevard #1'
  b_spb <- data$'Boulevard #2'
  p_msk <- data$'Park #1'
  p_spb <- data$'Park #2'
  
  six.colors <- c("Highway #1" = "#D49D6A", "Highway #2" = "#804815",
                  "Boulevard #1" = "#427A82", "Boulevard #2" = "#0E464E", 
                  "Park #1" = "#90BF60", "Park #2" = "#437313")
  
  Nsubs <- length(h_msk)
  df <- data.frame('Qscore'=c(h_msk,h_spb,b_msk,b_spb,p_msk,p_spb),
                   'cond'=c(rep('Highway #1', Nsubs),
                            rep('Highway #2', Nsubs),
                            rep('Boulevard #1', Nsubs),
                            rep('Boulevard #2', Nsubs),
                            rep('Park #1', Nsubs),
                            rep('Park #2', Nsubs)),
                   'id'=c(rep(1:Nsubs,6)))
  df$cond <- factor(df$cond,levels=c('Highway #1','Highway #2','Boulevard #1','Boulevard #2','Park #1','Park #2'))
  
  
  # Perform statistical tests
  
  ## Paired Wilcoxon & Kruskal-Wallis tests ####
  
  # We investigate the effect of the within-subject factor "Cond"
  # on the Questionnaire score
  
  # Kruskal-Wallis test 1
  # one.way.np <- df %>%
  #   kruskal_test(Qscore ~ cond)
  # one.way.np
  kruskal_test_result <- df %>%
    kruskal_test(Qscore ~ cond)
  print(kruskal_test_result)
  
  # Effect size of the Cond on Qscore
  # df %>%
  #   kruskal_effsize(Qscore ~ cond)
  effect_size_result <- df %>%
    kruskal_effsize(Qscore ~ cond)
  print(effect_size_result)
  
  # Wilcoxon test
  Wlx.test <- df %>% 
    wilcox_test(Qscore ~ cond, paired = TRUE,
                p.adjust.method = "BH", detailed = TRUE) %>%
    add_significance()
  print(Wlx.test)
  print(Wlx.test$conf.high)

  
  # Effect size for the pairwise comparisons
  # df %>% wilcox_effsize(Qscore ~ cond)
  pairwise_wilcox_effect_size <- df %>%
    wilcox_effsize(Qscore ~ cond)
  print(pairwise_wilcox_effect_size, n = 15, width = Inf)


  # For long plot titles not to be cut we need a function 
  # that inserts a newline character every certain number of characters
  wrap_title <- function(title) {
    paste(strwrap(title, width = 48), collapse = "\n")
  }

  # New version of the plotting function — plots mean scores with SEM (Standard Error of the Mean)
  # Calculate the mean and standard error for each condition
  df_summary <- df %>%
    group_by(cond) %>%
    summarise(mean_Qscore = mean(Qscore), sem = sd(Qscore) / sqrt(n()), .groups = 'drop')
  
  # Create a bar plot with ggplot2
  plot_env <- ggplot(df_summary, aes(x = cond, y = mean_Qscore, fill = cond)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    geom_errorbar(aes(ymin = mean_Qscore - sem, ymax = mean_Qscore + sem), width = 0.2) +
    geom_text(aes(label = round(mean_Qscore, 1), y = mean_Qscore + sem), vjust = -0.5, nudge_y = 0.3, size = 7.5) +
    scale_fill_manual(values = six.colors) +
    labs(
          title = wrap_title(plot_title),
          y = 'Mean score',
          #subtitle = "Comparison of questionnaires responses on a 1-7 scale",
          subtitle = "",
          #caption = paste("Kruskal-Wallis test: ", "p = ", one.way.np$p, ". pwc: Wilcoxon test")
          caption = "") +
        theme_minimal() + # White background
        theme(text = element_text(size = 15, family = "Arial"),
              plot.title = element_text(size = 48), # Set the title size here
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 24, hjust = 0.5, color = "black"),
              axis.title.y = element_text(size = 24, vjust = 3),
              axis.text.y = element_text(size = 10),
              plot.margin = unit(c(36, 36, 36, 36), "pt"),
              panel.grid.major = element_blank(), # Remove major grid lines
              panel.grid.minor = element_blank(), # Remove minor grid lines
              panel.border = element_blank(), # Remove panel border
              axis.line = element_line(color = "black"), # Add x and y axes
              legend.position = "none") +
        ylim(0, 7.5)
  
  # The previous version of the plotting function — plots median scores without error bars
  # Create a bar plot
  # plot_env <- ggbarplot(df, x = "cond", y = "Qscore", ylab = 'Median score',
  #     fill = "cond",
  #     error.plot = 'upper_errorbar',
  #     position = position_dodge(width = 0.8),
  #     add = c("median"), add.params = list(size = 1)) +
  #     geom_text(aes(label = round(..y.., 1)), stat = 'summary', fun = median, 
  #               nudge_y = 0.3, size = 7.5, vjust = -0.5) +
  #     #geom_jitter(size=1, width=0.1)+ # add subjects
  #     # add pairwise comparisons p-values
  #     #stat_pvalue_manual(Wlx.test, y.position = 10,
  #     #tip.length = 0, size=10,
  #     #step.increase = 0.2, hide.ns = T)+
  #     scale_fill_manual(values = six.colors) +
  #     # add labels
  #     labs(
  #       title = wrap_title(plot_title),
  #       #subtitle = "Comparison of questionnaires responses on a 1-7 scale",
  #       subtitle = "",
  #       #caption = paste("Kruskal-Wallis test: ", "p = ", one.way.np$p, ". pwc: Wilcoxon test")
  #       caption = "") +
  #     theme(text = element_text(size = 15, family = "Arial"),
  #           plot.title = element_text(size = 48), # Set the title size here
  #           axis.title.x = element_blank(),
  #           axis.text.x = element_text(size = 24, hjust = 0.5),
  #           axis.title.y = element_text(size = 24, vjust = 3),
  #           axis.text.y = element_text(size = 10),
  #           plot.margin = unit(c(36, 36, 36, 36), "pt"),
  #           theme(legend.position = "none")) +
  #     ylim(0, 7.5)
  
  # Show plot
  print(plot_env)
  
  return(plot_env)
}

# Array of corresponding plot titles
plot_titles <- c("“Poor Ecology of This Environment Threatens Residents' Health”", 
                 "“This Environment Is Very Boring”", 
                 "“It Is Impossible to Relax in This Environment”",
                 "“Do You Like This Environment?”",
                 "“This Environment Looks Pleasant”",
                 "“It's Interesting to Observe This Environment”",
                 "“This Environment Has Enough Green Spaces”",
                 "“Green Spaces in This Environment Are Relaxing”",
                 "“I Feel Safe in This Environment”",
                 "“This Environment Is Too Crowded”",
                 "“This Environment Has a Calm Atmosphere”",
                 "“Living in This Environment Is Quite Tiring”",
                 "“There Is Plenty to Do in This Environment”",
                 "“This Environment Is Too Noisy”",
                 "“Car Traffic in This Environment Is Very Annoying”",
                 "“This Environment Is Quiet”",
                 "“This Environment Is Clean”")

# Array of file names for saving
file_names_env <- c("1.poor_ecology", "2.very_boring", "3.impossible_relax", "4.like_environment", "5.looks_pleasant", "6.observe_environment", "7.green_spaces", "8.relaxing_spaces", "9.safe_environment", "10.crowded_environment", "11.calm_atmosphere", "12.tiring_living", "13.plenty_to_do", "14.noisy_environment", "15.car_traffic", "16.quiet_environment", "17.clean_environment")

# Construct the full file path
file_path_env <- "<your-path-to-data>"

# Array of sheet names or indices
sheets <- excel_sheets(file_path_env)

# Loop through the file names and corresponding plot titles
for (i in seq_along(sheets)) {
  
  # Run the function with the current file path and corresponding plot title
  current_plot <- questionnaire_environment_plot_xlsx(file_path_env, sheets[i], plot_titles[i])
  
  # Save the plot with the corresponding name
  output_directory_env <- "<your-path-to-data>"
  ggsave(paste0(output_directory_env, file_names_env[i], "_4800x3000.jpg"), current_plot, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")
}



# PLOTTING QUESTIONNAIRE GENERAL ANSWERS ----

# Array of file names
file_names_gen <- c("green_spaces", "quiet_pace", "noise_absence", "driving")
file_names_gen2 <- c("architecture", "cycling", "green_zones", "quiet_zones", "nature_contact", "water_proximity", "public_spaces", "urban_improvement", "local_community", "diversity", "safety", "school_access", "medical_access", "sports", "cultural_life", "supermarkets", "shopping_centers", "public_transport", "fast_paced_life", "ecology", "traffic_absence", "cleanliness", "environment_care", "emotional_connection")

# Function to process and plot data from an .xlsx file
questionnaire_general_plot_xlsx <- function(file_path, output_directory, output_directory2, file_names, file_names2) {
  
  # Read excel file
  data <- read_excel(file_path)
  
  # Loop through each column and plot
  for (i in seq_along(names(data))) {
    
    # Get column data
    column_data <- data[[i]]
    
    # Compute frequency of each response
    freq_data <- data.frame(table(factor(column_data, levels = 1:7)))
    names(freq_data) <- c('Response', 'Frequency')
    
    # Ensure all levels from 1 to 7 are present
    freq_data <- freq_data %>% 
      mutate(Response = as.integer(as.character(Response))) %>% 
      complete(Response = 1:7, fill = list(Frequency = 0))
    
    # For long plot titles not to be cut we need a function 
    # that inserts a newline character every certain number of characters
    wrap_title <- function(title) {
      paste(strwrap(title, width = 48), collapse = "\n")
    }
    
    # Plotting
    plot <- ggplot(freq_data, aes(x = factor(Response, levels = 1:7), y = Frequency)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = Frequency), vjust = -0.3, color = "black", size = 7.5) +
      labs(title = wrap_title(names(data)[i]),
           x = "Response scale from 1 to 7", 
           y = "Number of answers (30 in total)") +
      theme_minimal() +
      #ylim(0, 21) +
      scale_y_continuous(breaks = c(10, 20), limits = c(0, 21)) +
      theme(text = element_text(size = 15, family = "Arial"),
            plot.title = element_text(size = 48), 
            axis.title.x = element_text(size = 24, vjust = -3),
            axis.text.x = element_text(size = 24, hjust = 0.5),
            axis.title.y = element_text(size = 24, vjust = 3),
            axis.text.y = element_text(size = 24),
            plot.margin = unit(c(36, 36, 36, 36), "pt"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    
    # Show plot
    print(plot)
    
    # Save plot to the appropriate directory with the appropriate file name
    if (i <= 4) {
      ggsave(paste0(output_directory, file_names[i], "_4800x3000.jpg"), 
             plot, 
             width = 1280 / 80, 
             height = 800 / 80, 
             units = "in", 
             dpi = 300, 
             device = "jpeg")
    } else {
      ggsave(paste0(output_directory2, file_names2[i-4], "_4800x3000.jpg"), 
             plot, 
             width = 1280 / 80, 
             height = 800 / 80, 
             units = "in", 
             dpi = 300, 
             device = "jpeg")
    }
  }
}

# Construct the full file path
file_path_general <- "<your-path-to-data>"
output_directory_general <- "<your-path-to-data>"
output_directory_general_other <- "<your-path-to-data>"

# Run the function with the current file path
questionnaire_general_plot_xlsx(file_path_general, output_directory_general, output_directory_general_other, file_names_gen, file_names_gen2)
