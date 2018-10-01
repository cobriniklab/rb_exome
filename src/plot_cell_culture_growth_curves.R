#!/usr/bin/Rscript


# load required libraries -------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(grid)
library(cowplot)


# load required functions -------------------------------------------------

plot_cell_lines <- function(exome_curves, cohort, status){
  # browser()
  breaks = c(1, 3, 10, 30, 100)
  # plot_group <- filter(exome_curves, culture_group == culture_group)
  # plot_group <- plot_group[complete.cases(plot_group),]
  
  plot_group <- filter(exome_curves, culture_group == cohort, exome_status == status) %>% 
    group_by(sample) %>% 
    mutate(max = max(no_of_wells, na.rm = TRUE)) %>% 
    filter(max > 3)
  
  cl_plot <- ggplot(plot_group, aes(x=days_in_culture, y=no_of_wells, color = sample)) + 
    geom_line() + scale_y_log10(breaks = breaks, labels = breaks) + 
    theme(text = element_text(size=10), plot.title = element_text(lineheight=.5)) + 
    xlab("Days in Culture") + ylab("Number of Wells")
  
  cl_plot_facet <- cl_plot + facet_wrap(~sample)
  
  print(cl_plot)
  print(cl_plot_facet) 
  
  plot_path = paste0("~/rb_pipeline/doc/rb_cell_line_growth_curves/", cohort, "_", status, "_plot.png")
  plot_facet_path = paste0("~/rb_pipeline/doc/rb_cell_line_growth_curves/", cohort, "_", status, "_plot_facet.png")
  
  print(paste0("saving plots to: ", plot_path, " and ", plot_facet_path))
  
  save_plot(plot_facet_path, cl_plot_facet, base_aspect_ratio = 3.0)
  save_plot(plot_path, cl_plot, base_aspect_ratio = 1.5)
}

# load data ---------------------------------------------------------------

growth_curves <- read.table("~/rb_pipeline/data/rb_cell_line_growth_curves/growth_curves_20160125_update_20170825.csv", header = TRUE, sep = ",")

keep_cells <- scan("~/rb_pipeline/data/rb_cell_line_growth_curves/exomed_cell_lines.txt", what = character())
keep_cells <- gsub("-.*$", "", keep_cells)

growth_curves <- rownames_to_column(growth_curves, "days_in_culture") %>% 
  gather("sample", "no_of_wells", -days_in_culture) %>% 
  mutate(days_in_culture = as.integer(days_in_culture)) %>% 
  mutate(no_of_wells = as.integer(no_of_wells)) %>% 
  mutate(culture_group = ifelse(grepl("VC", sample), "vision_center", "reynolds")) %>%
  dplyr::mutate(cell_id = gsub("^.*\\.", "", sample)) %>% 
  dplyr::mutate(exome_status = ifelse((cell_id %in% keep_cells), "exome", "non-exome")) %>% 
  identity()


# plot data ---------------------------------------------------------------

plot_cell_lines(growth_curves, "vision_center", "exome")
plot_cell_lines(growth_curves, "reynolds", "exome")
plot_cell_lines(growth_curves, "vision_center", "non-exome")
plot_cell_lines(growth_curves, "reynolds", "non-exome")


# breaks = c(1, 3, 10, 30, 100)
# 
# # plot cell lines which were submitted for exome sequencing ---------------
# 
# vision_center <- filter(exomed_curves, culture_group == "vision_center")
# vision_center <- vision_center[complete.cases(vision_center),]
# 
# reynolds <- filter(exomed_curves, culture_group == "reynolds") %>% 
#   group_by(sample) %>% 
#   mutate(max = max(no_of_wells, na.rm = TRUE)) %>% 
#   filter(max > 3)
# 
# vision_center <- filter(exomed_curves[exomed_curves$culture_group == "vision_center",])
# 
# vc_plot <- ggplot(vision_center, aes(x=days_in_culture, y=no_of_wells, color = sample)) + 
#   geom_line() + scale_y_log10(breaks = breaks, labels = breaks) + 
#   theme(text = element_text(size=10), plot.title = element_text(lineheight=.5)) + 
#   xlab("Days in Culture") + ylab("Number of Wells")
# 
# vc_plot_facet <- vc_plot + facet_wrap(~sample)
# 
# print(vc_plot)
# print(vc_plot_facet)
# 
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/vc_exome_plot_facet.png", vc_plot_facet, base_aspect_ratio = 3.0)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/vc_exome_plot.png", vc_plot, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/doc/rb_cell_line_growth_curves/vc_exome_plot_facet.png", vc_plot_facet, base_aspect_ratio = 3.0)
# save_plot("~/rb_pipeline/doc/rb_cell_line_growth_curves/vc_exome_plot.png", vc_plot, base_aspect_ratio = 1.5)
# 
# reynolds_plot <- ggplot(reynolds, aes(x=days_in_culture, y=no_of_wells, color = sample)) + 
#   geom_line() +  scale_y_log10(breaks = breaks, labels = breaks) + 
#   theme(text = element_text(size=10), plot.title = element_text(lineheight=.5)) + 
#   xlab("Days in Culture") + ylab("Number of Wells")
# 
# print(reynolds_plot)
# reynolds_plot_facet <- reynolds_plot + facet_wrap(~sample)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/reynolds_exome_plot.png", reynolds_plot, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/reynolds_exome_plot_facet.png", reynolds_plot_facet, base_aspect_ratio = 3.0)

combined_plot <- plot_grid(vc_plot, reynolds_plot, align = "hv", nrow =2) + guides(colour = guide_legend(nrow = 2))
print(combined_plot)
ggsave("~/rb_pipeline/results/rb_cell_line_growth_curves/combined_exome_plot_20170914.pdf", combined_plot, width = 6, height = 8)


# # plot cell lines that were not submitted for exome sequencing ---------------

# vision_center <- filter(non_exomed_curves, culture_group == "vision_center")
# 
# reynolds <- filter(non_exomed_curves, culture_group == "reynolds") %>% 
#   group_by(sample) %>% 
#   mutate(max = max(no_of_wells, na.rm = TRUE)) %>% 
#   filter(max > 3)
# 
# vc_plot <- ggplot(vision_center, aes(x=days_in_culture, y=no_of_wells, color = sample)) + 
#   geom_line() + scale_y_log10(breaks = breaks, labels = breaks) + 
#   theme(text = element_text(size=10), plot.title = element_text(lineheight=.5)) + 
#   ggtitle("Cell Lines established under Vision Center") 
# 
# vc_plot_facet <- vc_plot + facet_wrap(~sample)
# 
# print(vc_plot)
# print(vc_plot_facet)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/vc_non_exome_plot.png", vc_plot, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/doc/rb_cell_line_growth_curves/vc_non_exome_plot.png", vc_plot, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/vc_non_exome_plot_facet.png", vc_plot_facet, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/doc/rb_cell_line_growth_curves/vc_non_exome_plot_facet.png", vc_plot_facet, base_aspect_ratio = 1.5)
# 
# reynolds_plot <- ggplot(reynolds, aes(x=days_in_culture, y=no_of_wells, color = sample)) + 
#   geom_line() +  scale_y_log10(breaks = breaks, labels = breaks) + 
#   theme(text = element_text(size=10), plot.title = element_text(lineheight=.5)) +
#   ggtitle("Cell Lines established under Reynolds Lab")
# 
# reynolds_plot_facet <- reynolds_plot + facet_wrap(~sample)
# 
# print(reynolds_plot)
# print(reynolds_plot_facet)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/reynolds_non_exome_plot.png", reynolds_plot, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/doc/rb_cell_line_growth_curves/reynolds_non_exome_plot.png", reynolds_plot, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/results/rb_cell_line_growth_curves/reynolds_non_exome_plot_facet.png", reynolds_plot_facet, base_aspect_ratio = 1.5)
# save_plot("~/rb_pipeline/doc/rb_cell_line_growth_curves/reynolds_non_exome_plot_facet.png", reynolds_plot_facet, base_aspect_ratio = 1.5)

combined_plot <- plot_grid(vc_plot, reynolds_plot, align = "hv", nrow =2) + guides(colour = guide_legend(nrow = 2))
print(combined_plot)
ggsave("~/rb_pipeline/results/rb_cell_line_growth_curves/combined_non_exome_plot_20170914.pdf", combined_plot, width = 6, height = 8)
