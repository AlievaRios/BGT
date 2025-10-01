library(dplyr)
library(ggplot2)
library(yaml)
library(optparse)

### Set to TRUE if you want to run import and processing even if file already exists
force_redo=TRUE
tracks_provided=NULL

BGT_dir <- NULL

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BGT_config file here if running the code in RStudio !!!!!!
  ### Demo path
  BGT_dir = paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/")
  pars = yaml.load_file("config_template.yml")
} else {
  ### Define a default BGT_dir for non-interactive sessions if not set
  ### You may adjust this path as necessary for your non-interactive environment
  if(is.null(BGT_dir)) {
    BGT_dir <- "~/BGT/"  # Adjust this path as necessary
  }
  
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BGT config file", metavar="character"),
    make_option(c("-f", "--force_redo"), action="store_true", default=FALSE, 
                help="Force the pipeline to re-import data even if files exists"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, 
                help="(Optional) Path to RDS file containing processed T cell track data", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  pars = yaml.load_file(opt$config)
  force_redo=opt$force_redo
  tracks_provided=opt$tracks_rds
}

### Setting data directory (if specified) and creating output directories
output_dir=paste0(BGT_dir, "/Results/Module1_Evaluation_of_Super_engager_population_dynamics_in_co_culture/")
dir.create(output_dir, recursive=TRUE)

imaging_time = pars$imaging_time
# Load data
master_clust_Live <- readRDS(pars$classified_tcell_track_data_filepath_rds)

# Data manipulations
# Group data by T cell line, cluster, and time, then calculate the count of records in each group
# Subsequently, regroups by T cell line and time to calculate the percentage of each group over the total count per timepoint.
Raw_per_overall <- master_clust_Live %>%
  group_by(tcell_line, cluster, Time) %>%
  summarise(n = n()) %>%
  filter(n() > 5) %>%
  group_by(tcell_line, Time) %>%
  mutate(perc = n * 100 / sum(n))

# Time conversion
Raw_per_overall$Time <- Raw_per_overall$Time / 30

# Filter out groups with fewer than 5 cells
#Raw_per_overall <- Raw_per_overall %>% filter(n > 5)

# Filter for imaging time, remove all the points before imaging time
Raw_per_overall <- Raw_per_overall %>%
  filter(Time >= imaging_time)

### Subset cluster 9
Raw_per_overall <- subset(Raw_per_overall, cluster == "9")

# Get unique T cell line names
unique_tcells <- unique(Raw_per_overall$tcell_line)
print(unique_tcells)

# Further subset the data for specific T cell lines
Raw_per_overall <- subset(Raw_per_overall, 
                          tcell_line %in% unique(Raw_per_overall$tcell_line))

imaging_time <- min(Raw_per_overall$Time)


# window offsets in Time units
offset_left  <- 15/30
offset_right <- 30/30

# palette keyed to tcell_line
lvls <- sort(unique(Raw_per_overall$tcell_line))
pal  <- brewer.pal(max(3, length(lvls)), "Set2")[seq_along(lvls)]
names(pal) <- lvls

library(scales)

x_min <- floor(min(Raw_per_overall$Time, na.rm = TRUE))
x_max <- ceiling(max(Raw_per_overall$Time, na.rm = TRUE))

# one window per tcell_line:
# take the Time at max(perc) within each line, then average if duplicates exist
win_df <- Raw_per_overall %>%
  group_by(tcell_line) %>%
  summarise(center_time = mean(Time[which(perc == max(perc, na.rm = TRUE))], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(window_start = center_time - offset_left,
         window_end   = center_time + offset_right)

Per2 <- ggplot(Raw_per_overall,
               aes(Time, perc, group = tcell_line, color = tcell_line)) +
  geom_smooth(size = 1, span = 0.5) +
  
  # shaded window per line (fill mapped to tcell_line so legend matches; hide fill legend)
  geom_rect(data = win_df,
            aes(xmin = window_start, xmax = window_end,
                ymin = -Inf, ymax = Inf, fill = tcell_line),
            alpha = 0.1, inherit.aes = FALSE, show.legend = FALSE) +
  
  # verticals in the SAME color as each line; no extra legend keys
  geom_vline(data = win_df, aes(xintercept = window_start, color = tcell_line),
             linetype = "dotdash", show.legend = FALSE) +
  geom_vline(data = win_df, aes(xintercept = window_end,   color = tcell_line),
             linetype = "dotdash", show.legend = FALSE) +
  
  theme_bw() +
  ylab("% T cells in cluster 9 (super-engagers)") +
  xlab("Time in Co-culture (Hours)") +
  scale_x_continuous(breaks = c(imaging_time,
                                seq(from = ceiling(imaging_time),
                                    to = max(Raw_per_overall$Time, na.rm = TRUE), by = 1))) +
  scale_color_manual(values = pal, name = "T Cell Type") +
  scale_fill_manual(values  = pal, guide = "none") +
  ggtitle("T-cell engagement at different time points")+
  scale_x_continuous(
    breaks = seq(x_min, x_max, by = 1),                # integers
    labels = label_number(accuracy = 1),               # show as 0,1,2,...
    minor_breaks = seq(x_min, x_max, by = 0.25)        # quarters
  )

ggsave(filename = paste0(output_dir, "T-cell_engagementVsTime.png"), plot = Per2, width = 8, height = 6)

# Save data
write.csv(Raw_per_overall, file.path(paste0(output_dir, "T-cell_engagementVsTime.csv")), row.names = FALSE)


# Display plots
Per2

