# Load some packages that we will need
library(rstatix)
library(ggplot2)
library(ggpubr)
library(R.matlab)

library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(dplyr)

library(effectsize)


### Data Loading & Preparation ====

## Set the working directory with the data for statistical analysis
setwd("<your-path-to-data>")

## Set the output directory to save plots
output_directory <- "<your-path-to-data>"
output_directory_freq <- "<your-path-to-data>"


## Function to load data based on the keyword
load_data <- function(keyword) {
  
  # Video names
  videos <- c('h_msk', 'h_spb', 'b_msk', 'b_spb', 'p_msk', 'p_spb', 'sea', 'mov')
  
  # Initialize list to store the data
  data_list <- list()
  
  # Loop over each video name
  for (i in seq_along(videos)) {
    # Construct the file name
    file_name <- paste0(videos[i], '_', keyword, '_ISC_persubject_sum.mat')
    
    # Load the data
    data_list[[i]] <- as.vector(readMat(file_name)$ISC.persubject)
  }
  
  return(data_list)
}


## GENERAL vs. FREQUENCY BANDS — pass the desired input to the function:
## "all" — for the general case
## "delta"/"theta"/"alpha"/"beta"/"gamma" — for frequency bands 
data <- load_data('all')


## Assign data to variables
h_msk <- data[[1]]
h_spb <- data[[2]]
b_msk <- data[[3]]
b_spb <- data[[4]]
p_msk <- data[[5]]
p_spb <- data[[6]]
sea <- data[[7]]
mov <- data[[8]]


## set fixed colors for the environments from the palette https://paletton.com/#uid=73l210kllllaFw0g0qFqFg0w0aF
three.colors <- c("Highway" = "#AA6F39", "Boulevard" = "#246068", "Park" = "#669933")
sm.colors <- c("Control #1" = "#BA5D7E", "Control #2" = "#953255")
three_sm.colors <- c("Control #1" = "#BA5D7E","Highway" = "#AA6F39", "Boulevard" = "#246068", "Park" = "#669933","Control #2" = "#953255")
three_s.colors <- c("Control #1" = "#BA5D7E","Highway" = "#AA6F39", "Boulevard" = "#246068", "Park" = "#669933")
msk.colors <- c("Highway #1" = "#D49D6A", "Boulevard #1" = "#427A82","Park #1" = "#90BF60")
spb.colors <- c("Highway #2" = "#804815", "Boulevard #2" = "#0E464E", "Park #2" = "#437313")
six.colors <- c("Highway #1" = "#D49D6A", "Boulevard #1" = "#427A82",
                "Highway #2" = "#804815", "Boulevard #2" = "#0E464E", 
                "Park #2" = "#437313", "Park #1" = "#90BF60")
all.colors <- c("Control #1" = "#BA5D7E", "Highway #1" = "#D49D6A", "Boulevard #1" = "#427A82",
                "Highway #2" = "#804815", "Boulevard #2" = "#0E464E", 
                "Park #2" = "#437313", "Park #1" = "#90BF60", "Control #2" = "#953255")


### Since we don't have different groups with different stimuli
### we perform only the WITHIN-SUBJECT ANALYSIS - 
### for a within-subject factor Condition (video)


### FUNCTIONS ----

## Function to check statistical assumptions
check_assumptions <- function(df, within_factors) {
  
  # Check outliers 
  df %>%
    group_by(across(all_of(within_factors))) %>%
    rstatix::identify_outliers(ISC) # should be no extreme outliers
  
  # Check normality - data should be normally distributed
  # Shapiro-Wilk test
  print(df %>%
    group_by(across(all_of(within_factors))) %>%
    rstatix::shapiro_test(ISC)) #p should be >0.05
  
  # Draw a Q-Q plot by cond
  plot <- ggpubr::ggqqplot(df, x = "ISC", facet.by = within_factors)
  print(plot)
  # the Q-Q plot should follow the line closely 
  
  # Check the equality of variances — Levene's Test for Homogeneity of Variance
  # In the context of an ANOVA and t-test, variances of our groups 
  # should not significantly different from each other (p > 0.05)
  formula <- as.formula(paste("ISC ~", paste(within_factors, collapse = "*")))
  df %>% rstatix::levene_test(formula) # p should be >0.05
}

## Function to prevent long titles from cutting by the plot_function.
## It inserts a newline character every certain number of characters -
wrap_title <- function(title) {
  paste(strwrap(title, width = 48), collapse = "\n")
}

## Plotting function - to make bar or box plots on demand
plot_function <- function(df, x, y, ylab_value, fill, fill_colors, 
                          pwc_data, y_position, title_text, subtitle_text, 
                          caption, x_text_size, hide_ns = T, 
                          add_jitter = FALSE, add_pvalue = TRUE, add_paired = FALSE,
                          y_limits = NULL) {
  # Plots a paired box plot if requested
  if(add_paired) {
    plot <- ggpaired(df, x = x, y = y, 
                     fill=fill_colors,
                     point.size=2)+
      ylab(ylab_value)
    
      if (!is.null(y_limits)) {
        plot <- plot + ylim(y_limits)  # Setting y limits only if y_limits is not NULL
      }
      
      plot <- plot +
      labs( # add labels
        title = wrap_title(title_text),
        subtitle = subtitle_text,
        caption = caption)+
      theme(text = element_text(size = 15, family = "Arial"),
            plot.title = element_text(size = 48), # Set the title size here
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = x_text_size, hjust = 0.5),
            axis.title.y = element_text(size = 24, vjust = 3),
            axis.text.y = element_text(size = 15),
            plot.margin = unit(c(36, 36, 36, 36), "pt"))
  } else {
    # Plots a bar plot by default
    plot <- ggbarplot(df, x = x, y = y, ylab = ylab_value,
                      fill = fill,
                      error.plot = 'upper_errorbar',
                      position = position_dodge(width = 0.8),
                      add = c("mean_se"), add.params = list(size = 1)) +
      scale_fill_manual(values = fill_colors)
      
      if (!is.null(y_limits)) {
        plot <- plot + ylim(y_limits)  # Setting y limits only if y_limits is not NULL
      }
      
      plot <- plot +
      # add labels
      labs(
        title = wrap_title(title_text),
        subtitle = subtitle_text,
        caption = caption) +
      theme(text = element_text(size = 15, family = "Arial"),
            plot.title = element_text(size = 48), # Set the title size here
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = x_text_size, hjust = 0.5),
            axis.title.y = element_text(size = 24, vjust = 3),
            axis.text.y = element_text(size = 15),
            plot.margin = unit(c(36, 36, 36, 36), "pt")) +
      theme(legend.position = "none")
  }
  
  # Adds pairwise comparisons p-values by default
  if(add_pvalue && !add_paired) {
    plot <- plot + stat_pvalue_manual(pwc_data, y.position = y_position,
                                      tip.length = 0, size = 10,
                                      step.increase = 0.1, hide.ns = hide_ns)
  }
  # Adds jitter (dots representing the subjects) if requested
  if(add_jitter && !add_paired) {
    plot <- plot + geom_jitter(size = 1, width = 0.1)
  }
  
  return(plot)
}

## Function to calculate statistics:
## - One- or Two-way repeated measures ANOVA
## - T-test or Wilcoxon test — for pairwise comparisons

calculate_statistics <- function(df, within_factors, t_test_factor, label_position_factor, paired = TRUE, wilcox = FALSE, conf.level = 0.95) {
  
  # Create the formula for the tests dynamically based on the factors
  formula <- as.formula(paste("ISC ~", paste(t_test_factor)))
  
  # One-way or two-way repeated measures Anova
  res.aov <- rstatix::anova_test(
    data = df, dv = ISC, wid = id,
    within = within_factors
  )
  anova_table <- rstatix::get_anova_table(res.aov)
  print(anova_table)
  
  # T-test - pairwise comparisons between conditions
  pwc <- df %>%
    rstatix::pairwise_t_test(
      formula, paired = paired,
      p.adjust.method = "bonferroni"
    )  
  print(pwc)
  
  # Wilcoxon test - if the assumptions for the data aren't met
  if (wilcox) {
    Wlx.test <- df %>%
      rstatix::wilcox_test(formula, paired = paired,
                           p.adjust.method = "bonferroni") %>%
      add_significance()
    print(Wlx.test)
    pwc <- Wlx.test
  }
  
  # Calculate the positions of text labels for significance levels 
  pwc <- pwc %>% add_xy_position(x = label_position_factor)
  
  list(anova_table = anova_table, pwc = pwc)
}

## Function to calculate statistics with Effect Sizes and Confidence Intervals:
calculate_statistics_w_ES_CI <- function(df, within_factors, t_test_factor, label_position_factor, paired = TRUE, wilcox = FALSE) {
  
  # Create the formula for the tests dynamically based on the factors
  formula <- as.formula(paste("ISC ~", paste(t_test_factor)))
  
  # One-way or two-way repeated measures Anova
  res.aov <- rstatix::anova_test(
    data = df, dv = ISC, wid = id,
    within = within_factors
  )
  anova_table <- rstatix::get_anova_table(res.aov)
  print(anova_table)
  
  # T-test - pairwise comparisons between conditions
  pwc <- df %>%
    rstatix::pairwise_t_test(
      formula, paired = paired,
      p.adjust.method = "bonferroni"
    )
  print(pwc)
  
  # Add confidence intervals and effect sizes for each pair
  pairs <- combn(levels(df$cond), 2, simplify = FALSE)
  for (pair in pairs) {
    data_pair <- df[df$cond %in% pair, ]
    t_test_result <- t.test(ISC ~ cond, data = data_pair, paired = paired)
    conf_int <- t_test_result$conf.int
    
    # Calculate Cohen's D — an effect size used to indicate 
    # the standardized difference between two means
    cohen_d_result <- effectsize::cohens_d(ISC ~ cond, data = data_pair, paired = paired)
    print(cohen_d_result)
    
    cohen_d <- cohen_d_result[1, "estimate"]
    cat("Pair:", pair, "Confidence interval:", conf_int, "Cohen's D:", cohen_d, "\n")
  }
  
  # Wilcoxon test - if the assumptions for the data aren't met
  if (wilcox) {
    Wlx.test <- df %>%
      rstatix::wilcox_test(formula, paired = paired, conf.int = TRUE,
                           p.adjust.method = "bonferroni") %>%
      add_significance()
    print(Wlx.test)
    pwc <- Wlx.test
  }
  
  # Calculate the positions of text labels for significance levels 
  pwc <- pwc %>% add_xy_position(x = label_position_factor)
  
  # Print the pwc data frame
  print(pwc)
  
  list(anova_table = anova_table, pwc = pwc)
}



### 6 URBAN VIDEOS ====

## 1. Create a data frame
Nsubs=length(h_msk)
df_six=data.frame('ISC'=c(h_msk,b_msk,h_spb,b_spb,p_spb,p_msk),
                     'cond'=c(rep('Highway #1', Nsubs),
                                   rep('Boulevard #1', Nsubs),
                                   rep('Highway #2', Nsubs),
                                   rep('Boulevard #2', Nsubs),
                                   rep('Park #2', Nsubs),
                                   rep('Park #1', Nsubs)),
                     'id'=c(rep(1:Nsubs,6)))
df_six$cond=factor(df_six$cond,levels=c('Highway #1','Boulevard #1','Highway #2','Boulevard #2','Park #2','Park #1'))


## 2.Check assumptions for Anova & T-tests

# Call the function
check_assumptions(df_six, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video)

# Call the plotting function to make the barplot with individual results
plot_dots_six <- plot_function(df_six, "cond", "ISC", "ISC", "cond", six.colors, pwc.six, 0.025, 
                          "Inter-Subject Neural Synchronization (ISC) Across Six Urban Environment Videos",
                          "", "", 24, add_jitter = TRUE, add_pvalue = FALSE,
                          y_limits = c(0, 0.04))
plot_dots_six
ggsave(paste0(output_directory, "Six_dots_4800x3000.jpg"), plot_dots_six, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the boxplot diagram with paired results
plot_paired_six <- plot_function(df_six, "cond", "ISC", "ISC", "cond", six.colors, pwc.six, 0.025, 
                          "Inter-Subject Neural Synchronization (ISC) Across Six Urban Environment Videos",
                          "", "", 24, add_paired = TRUE,
                          y_limits = c(0, 0.04))
plot_paired_six
ggsave(paste0(output_directory, "Six_paired_4800x3000.jpg"), plot_paired_six, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.six <- calculate_statistics2(df_six, c("cond"), "cond", "cond")
res.aov.six <- stats_results.six$anova_table
pwc.six <- stats_results.six$pwc

# To calculate statistics with Effect Sizes and Confidence Intervals:
#stats_results.six <- calculate_statistics_w_ES_CI(df_six, c("cond"), "cond", "cond")
#res.aov.six <- stats_results.six$anova_table
#pwc.six <- stats_results.six$pwc


## 5. Report & Visualization — the main bar plot with the results

# For the general case.
# Call the plotting function to create and save the first plot with short y label 'ISC'
plot_six <- plot_function(df_six, "cond", "ISC", "ISC", "cond", six.colors, pwc.six, 0.025, 
                           "Inter-Subject Neural Synchronization (ISC) Across Six Urban Environment Videos",
                           "", "", 24, y_limits = c(0, 0.04))
plot_six
ggsave(paste0(output_directory, "Six_anova_4800x3000.jpg"), plot_six, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")
#ggsave(paste0(output_directory, "Six_anova_4800x3000.svg"), plot_six, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "svg")


# For frequency bands — insert the required title:
# "Delta Frequency Band (0.5-4 Hz)"
# "Theta Frequency Band (4-8 Hz)"
# "Alpha Frequency Band (8-12 Hz)"
# "Beta Frequency Band (12-25 Hz)"
# "Gamma Frequency Band (25-55 Hz)"
plot_six_freq <- plot_function(df_six, "cond", "ISC", "ISC", "cond", six.colors, pwc.six, 0.06, 
                          "Theta Frequency Band (4-8 Hz)",
                          "", 
                          "", 24, y_limits = c(0, 0.105)) # set y_limits for the same scale
plot_six_freq
# Insert the required band name
ggsave(paste0(output_directory_freq, "Six_anova_", "theta", "_4800x3000.jpg"), plot_six_freq, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")




### SEA vs MOVIE - comparing control conditions ====
### to check whether we use the ISC method correctly

## 1. Create a data frame

df_sm = data.frame('ISC'=c(sea,mov),
                   'cond'=c(rep('Control #1', length(sea)),
                            rep('Control #2', length(mov))),
                   'id'=c(rep(1:length(mov),2)))


## 2.Check assumptions for a T-test

# Call the function to check statistical assumptions
check_assumptions(df_sm, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the barplot with individual results
plot_dots_sm <- plot_function(df_sm, "cond", "ISC", "ISC", "cond", sm.colors, Wlx.test.sm, 0.2, 
                               "Inter-Subject Neural Synchronization (ISC) Across Two Control Videos",
                               "", "", 24, add_jitter = TRUE, add_pvalue = FALSE)
plot_dots_sm
ggsave(paste0(output_directory, "Control_dots_4800x3000.jpg"), plot_dots_sm, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the boxplot diagram with paired results
plot_paired_sm <- plot_function(df_sm, "cond", "ISC", "ISC", "cond", sm.colors, Wlx.test.sm, 0.2, 
                                 "Inter-Subject Neural Synchronization (ISC) Across Two Control Videos",
                                 "", "", 24, add_paired = TRUE)
plot_paired_sm
ggsave(paste0(output_directory, "Control_paired_4800x3000.jpg"), plot_paired_sm, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# T-test for pairwise comparison between Movie and Sea
# Wilcoxon test - if the assumptions for the data aren't met
stats_results.sm <- calculate_statistics(df_sm, c("cond"), "cond", "cond", wilcox = TRUE)
Wlx.test.sm <- stats_results.sm$pwc


## 5. Report & Vusialization

# Call the plotting function to create and save the first plot with short y label 'ISC'
plot_sm <- plot_function(df_sm, "cond", "ISC", "ISC", "cond", sm.colors, Wlx.test.sm, 0.2, 
                            "Inter-Subject Neural Synchronization (ISC) Across Two Control Videos",
                            "", "", 24, y_limits = c(0, 0.2))
plot_sm
ggsave(paste0(output_directory, "Control_4800x3000.jpg"), plot_sm, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



### 8 VIDEOS: 6 urban + 2 control videos ====

## 1. Create a data frame

Nsubs=length(h_msk)
df_all=data.frame('ISC'=c(sea,h_msk,b_msk,h_spb,b_spb,p_spb,p_msk,mov),
                          'cond'=c(rep('Control #1', Nsubs),
                                        rep('Highway #1', Nsubs),
                                        rep('Boulevard #1', Nsubs),
                                        rep('Highway #2', Nsubs),
                                        rep('Boulevard #2', Nsubs),
                                        rep('Park #2', Nsubs),
                                        rep('Park #1', Nsubs),
                                        rep('Control #2', Nsubs)),
                          'id'=c(rep(1:Nsubs,8)))
df_all$cond=factor(df_all$cond,levels=c('Control #1','Highway #1','Boulevard #1','Highway #2','Boulevard #2','Park #2','Park #1','Control #2'))


## 2.Check assumptions for Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_all, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the bar plot with individual results
plot_dots_all <- plot_function(df_all, "cond", "ISC", "ISC", "cond", all.colors, pwc.all, 0.025, 
                               "Inter-Subject Neural Synchronization (ISC) Across Six Urban Environment and Two Control Videos",
                               "", "", 20, add_jitter = TRUE, add_pvalue = FALSE,
                               y_limits = c(0, 0.2))
plot_dots_all
ggsave(paste0(output_directory, "All_dots_4800x3000.jpg"), plot_dots_all, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the box plot diagram with paired results
plot_paired_all <- plot_function(df_all, "cond", "ISC", "ISC", "cond", all.colors, pwc.all, 0.025, 
                                 "Inter-Subject Neural Synchronization (ISC) Across Six Urban Environment and Two Control Videos",
                                 "", "", 24, add_paired = TRUE)
plot_paired_all
ggsave(paste0(output_directory, "All_paired_4800x3000.jpg"), plot_paired_all, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.all <- calculate_statistics(df_all, c("cond"), "cond", "cond")
res.aov.all <- stats_results.all$anova_table
pwc.all <- stats_results.all$pwc

# To calculate statistics with Effect Sizes and Confidence Intervals:
stats_results.all <- calculate_statistics_w_ES_CI(df_all, c("cond"), "cond", "cond")
res.aov.all <- stats_results.all$anova_table
pwc.all <- stats_results.all$pwc


## 5. Report & Visualization — the main bar plot with the results

# Call the plotting function to make the first plot with short ylab
plot_all <- plot_function(df_all, "cond", "ISC", "ISC", "cond", all.colors, pwc.all, 0.15, 
                          "Inter-Subject Neural Synchronization (ISC) Across Six Urban Environment and Two Control Videos",
                          "",
                          "", 20, y_limits = c(0, 0.35))
plot_all
ggsave(paste0(output_directory, "All_anova_4800x3000.jpg"), plot_all, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



### CITY #1 videos ====

## 1. Create a data frame
Nsubs=length(h_msk)
df_msk=data.frame('ISC'=c(h_msk,b_msk,p_msk),
                  'cond'=c(rep('Highway #1', Nsubs),
                                rep('Boulevard #1', Nsubs),
                                rep('Park #1', Nsubs)),
                  'id'=c(rep(1:Nsubs,3)))
df_msk$cond=factor(df_msk$cond,levels=c('Highway #1','Boulevard #1','Park #1'))


## 2.Check assumptions for Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_msk, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the bar plot with individual results
plot_dots_msk <- plot_function(df_msk, "cond", "ISC", "ISC", "cond", msk.colors, pwc.msk, 0.025, 
                               "Inter-Subject Neural Synchronization (ISC) Across the City #1 Urban Environment Videos",
                               "", "", 24, add_jitter = TRUE, add_pvalue = FALSE,
                               y_limits = c(0, 0.04))
plot_dots_msk
ggsave(paste0(output_directory, "Msk_dots_4800x3000.jpg"), plot_dots_msk, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the box plot diagram with paired results
plot_paired_msk <- plot_function(df_msk, "cond", "ISC", "ISC", "cond", msk.colors, pwc.msk, 0.025, 
                                 "Inter-Subject Neural Synchronization (ISC) Across the City #1 Urban Environment Videos",
                                 "", "", 24, add_paired = TRUE,
                                 y_limits = c(0, 0.04))
plot_paired_msk
ggsave(paste0(output_directory, "Msk_paired_4800x3000.jpg"), plot_paired_msk, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.msk <- calculate_statistics(df_msk, c("cond"), "cond", "cond")
res.aov.msk <- stats_results.msk$anova_table
pwc.msk <- stats_results.msk$pwc


## 5. Report & Visualization — the main bar plot with the results

# Call the plotting function to make the first plot with short ylab
plot_msk <- plot_function(df_msk, "cond", "ISC", "ISC", "cond", msk.colors, pwc.msk, 0.025, 
                           "Inter-Subject Neural Synchronization (ISC) Across the City #1 Urban Environment Videos",
                           "",
                          "", 24, y_limits = c(0, 0.04))
plot_msk
ggsave(paste0(output_directory, "Msk_anova_4800x3000.jpg"), plot_msk, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



### CITY #2 videos ====

## 1. Create data frame
Nsubs=length(h_spb)
df_spb=data.frame('ISC'=c(h_spb,b_spb,p_spb),
                  'cond'=c(rep('Highway #2', Nsubs),
                           rep('Boulevard #2', Nsubs),
                                rep('Park #2', Nsubs)),
                  'id'=c(rep(1:Nsubs,3)))
df_spb$cond=factor(df_spb$cond,levels=c('Highway #2','Boulevard #2','Park #2'))


## 2.Check assumptions for Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_spb, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the bar plot with individual results
plot_dots_spb <- plot_function(df_spb, "cond", "ISC", "ISC", "cond", spb.colors, pwc.spb, 0.025, 
                               "Inter-Subject Neural Synchronization (ISC) Across the City #2 Urban Environment Videos",
                               "", "", 24, add_jitter = TRUE, add_pvalue = FALSE,
                               y_limits = c(0, 0.04))
plot_dots_spb
ggsave(paste0(output_directory, "Spb_dots_4800x3000.jpg"), plot_dots_spb, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the box plot diagram with paired results
plot_paired_spb <- plot_function(df_spb, "cond", "ISC", "ISC", "cond", spb.colors, pwc.spb, 0.025, 
                                 "Inter-Subject Neural Synchronization (ISC) Across the City #2 Urban Environment Videos",
                                 "", "", 24, add_paired = TRUE,
                                 y_limits = c(0, 0.04))
plot_paired_spb
ggsave(paste0(output_directory, "Spb_paired_4800x3000.jpg"), plot_paired_spb, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.spb <- calculate_statistics(df_spb, c("cond"), "cond", "cond")
res.aov.spb <- stats_results.spb$anova_table
pwc.spb <- stats_results.spb$pwc


## 5. Report & Visualization — the main bar plot with the results

# Call the plotting function to make the first plot with short ylab
plot_spb <- plot_function(df_spb, "cond", "ISC", "ISC", "cond", spb.colors, pwc.spb, 0.025, 
                          "Inter-Subject Neural Synchronization (ISC) Across the City #2 Urban Environment Videos",
                          "",
                          "", 24, y_limits = c(0, 0.04))
plot_spb
ggsave(paste0(output_directory, "Spb_anova_4800x3000.jpg"), plot_spb, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



### 3 URBAN - NO CITY ====

### Here we collapse the data from the two cities for each urban environment

## 1. Create a data frame
## we take element-wise means for City #1 & City #2 data in each condition

Nsubs=length(h_msk)
df_three=data.frame('ISC'=c(rowMeans(cbind(h_msk, h_spb)),
                            rowMeans(cbind(b_msk, b_spb)),
                            rowMeans(cbind(p_msk, p_spb))),
                  'cond'=c(rep('Highway', Nsubs),
                           rep('Boulevard', Nsubs),
                           rep('Park', Nsubs)),
                  'id'=c(rep(1:Nsubs,3)))
df_three$cond=factor(df_three$cond,levels=c('Highway','Boulevard','Park'))


## 2.Check assumptions for Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_three, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the bar plot with individual results
plot_dots_three <- plot_function(df_three, "cond", "ISC", "ISC", "cond", three.colors, pwc.three, 0.025, 
                               "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities",
                               "", "", 24, add_jitter = TRUE, add_pvalue = FALSE,
                               y_limits = c(0, 0.04))
plot_dots_three
ggsave(paste0(output_directory, "Three_dots_4800x3000.jpg"), plot_dots_three, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the box plot diagram with paired results
plot_paired_three <- plot_function(df_three, "cond", "ISC", "ISC", "cond", three.colors, pwc.three, 0.025, 
                                 "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities",
                                 "", "", 24, add_paired = TRUE,
                                 y_limits = c(0, 0.04))
plot_paired_three
ggsave(paste0(output_directory, "Three_paired_4800x3000.jpg"), plot_paired_three, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.three <- calculate_statistics(df_three, c("cond"), "cond", "cond")
res.aov.three <- stats_results.three$anova_table
pwc.three <- stats_results.three$pwc


## 5. Report & Visualization — the main bar plot with the results

# Call the plotting function to make the first plot with short ylab
plot_three <- plot_function(df_three, "cond", "ISC", "ISC", "cond", three.colors, pwc.three, 0.025, 
                          "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities",
                          "",
                          "", 24, y_limits = c(0, 0.04))
plot_three
ggsave(paste0(output_directory, "Three_anova_4800x3000.jpg"), plot_three, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



#### 3 URBAN + 2 CONTROL videos ====

## 1. Create a data frame

df_three_sm=rbind(df_three,df_sm)
summary(df_three_sm)
df_three_sm$cond=factor(df_three_sm$cond,
                              levels=c('Control #1','Highway','Boulevard','Park','Control #2'))


## 2.Check assumptions for Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_three_sm, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the bar plot with individual results
plot_dots_three_sm <- plot_function(df_three_sm, "cond", "ISC", "ISC", "cond", three_sm.colors, pwc.three_sm, 0.025, 
                                 "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities and Two Control Videos",
                                 "", "", 24, add_jitter = TRUE, add_pvalue = FALSE,
                                 y_limits = c(0, 0.2))
plot_dots_three_sm
ggsave(paste0(output_directory, "Three_SM_dots_4800x3000.jpg"), plot_dots_three_sm, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the box plot diagram with paired results
plot_paired_three_sm <- plot_function(df_three_sm, "cond", "ISC", "ISC", "cond", three_sm.colors, pwc.three_sm, 0.025, 
                                   "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities and Two Control Videos",
                                   "", "", 24, add_paired = TRUE,
                                   y_limits = c(0, 0.2))
plot_paired_three_sm
ggsave(paste0(output_directory, "Three_SM_paired_4800x3000.jpg"), plot_paired_three_sm, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.three_sm <- calculate_statistics(df_three_sm, c("cond"), "cond", "cond")
res.aov.three_sm <- stats_results.three_sm$anova_table
pwc.three_sm <- stats_results.three_sm$pwc


## 5. Report & Visualization — the main bar plot with the results

# Call the plotting function to make the first plot with short ylab
plot_three_sm <- plot_function(df_three_sm, "cond", "ISC", "ISC", "cond", three_sm.colors, pwc.three_sm, 0.2, 
                            "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities and Two Control Videos",
                            "",
                            "", 24, y_limits = c(0, 0.35))
plot_three_sm
ggsave(paste0(output_directory, "Three_SM_anova_4800x3000.jpg"), plot_three_sm, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



### 3 URBAN + SEA videos ====

## 1. Create a data frame

df_s=data.frame('ISC'=c(sea),
                'cond'=c(rep('Control #1', length(sea))),
                'id'=c(rep(1:length(sea),1)))

df_three_s=rbind(df_three,df_s)
summary(df_three_s)
df_three_s$cond=factor(df_three_s$cond,
                              levels=c('Control #1','Highway',
                                       'Boulevard','Park'))


## 2.Check assumptions for Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_three_s, c("cond")) # p should be >0.05


## 3. Visualize the distribution of data in each Condition (video) 

# Call the plotting function to make the bar plot with individual results
plot_dots_three_s <- plot_function(df_three_s, "cond", "ISC", "ISC", "cond", three_s.colors, pwc.three_s, 0.025, 
                                    "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities and One Control Video",
                                    "", "", 24, add_jitter = TRUE, add_pvalue = FALSE,
                                   y_limits = c(0, 0.04))
plot_dots_three_s
ggsave(paste0(output_directory, "Three_S_dots_4800x3000.jpg"), plot_dots_three_s, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")

# Call the plotting function to make the box plot diagram with paired results
plot_paired_three_s <- plot_function(df_three_s, "cond", "ISC", "ISC", "cond", three_s.colors, pwc.three_s, 0.025, 
                                      "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities and One Control Video",
                                      "", "", 24, add_paired = TRUE,
                                     y_limits = c(0, 0.04))
plot_paired_three_s
ggsave(paste0(output_directory, "Three_S_paired_4800x3000.jpg"), plot_paired_three_s, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")


## 4. Calculate statistics

# Call the function to calculate:
# - One-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.three_s <- calculate_statistics(df_three_s, c("cond"), "cond", "cond")
res.aov.three_s <- stats_results.three_s$anova_table
pwc.three_s <- stats_results.three_s$pwc


## 5. Report & Visualization — the main bar plot with the results

# Call the plotting function to make the first plot with short ylab
plot_three_s <- plot_function(df_three_s, "cond", "ISC", "ISC", "cond", three_s.colors, pwc.three_s, 0.025, 
                               "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos Collapsed by Cities and One Control Video",
                               "",
                               "", 24, y_limits = c(0, 0.04))
plot_three_s
ggsave(paste0(output_directory, "Three_S_anova_4800x3000.jpg"), plot_three_s, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")



### VIDEOS & CITIES - two-way Anova ====
### Here we perform a within-subjects analysis of two conditions: video and city

## 1. Create data frame
Nsubs=length(h_msk)
df_city=data.frame('ISC'=c(h_msk,b_msk,p_msk,h_spb,b_spb,p_spb),
                   'video'=c(rep(c(rep('Highway', Nsubs),
                                   rep('Boulevard', Nsubs),
                                   rep('Park', Nsubs)),2)),
                   'city'=c(rep('City #1',Nsubs*3),rep('City #2',Nsubs*3)),
                   'id'=c(rep(1:Nsubs,6)))
df_city$video=factor(df_city$video,levels=c('Highway','Boulevard','Park'))


## 2.Check assumptions for two-way Anova & T-tests

# Call the function to check statistical assumptions
check_assumptions(df_city, c("video", "city")) # p should be >0.05


## 3. Calculate statistics

# Call the function to calculate:
# - Two-way repeated measures Anova
# - T-test for pairwise comparisons between conditions
stats_results.city <- calculate_statistics(df_city, c("video", "city"), "video", "city")
res.aov.city <- stats_results.city$anova_table
pwc.city <- stats_results.city$pwc


## 4. Report & Visualization — the main bar plot with the results

# Call the plotting function to create and save the first plot with short y label 'ISC'
plot_city <- plot_function(df_city, "city", "ISC", "ISC", "video", three.colors, pwc.city, 0.025, 
                          "Inter-Subject Neural Synchronization (ISC) Across Urban Environment Videos in Two Cities",
                          "", 
                          "", 24,
                          y_limits = c(0, 0.04))
plot_city
ggsave(paste0(output_directory, "City_2way_anova_4800x3000.jpg"), plot_city, width = 1280/80, height = 800/80, units = "in", dpi = 300, device = "jpeg")
