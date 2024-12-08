library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidydr)
library(extrafont)
library(ggplot2)

font_import(prompt = F)
loadfonts(device = "pdf") 
loadfonts(device = "postscript") 

my_plot_theme <- function(legend = "right", base_family = "Arial", x.text.angle = 0) {
  theme_pubr(legend = legend, base_family = base_family, x.text.angle = x.text.angle) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme(
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      plot.title = element_text(size = 8, face = "bold", margin = margin(t = 0)),
      plot.subtitle = element_text(size = 7, margin = margin(t = 0)),
      plot.caption = element_text(size = 6),
      plot.margin = margin(2, 2, 2, 2),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent"),
    )
}

my_sc_plot_theme <- function(xlength = 0.25, ylength = 0.25,
                             arrow = grid::arrow(length = unit(0.05, "inches"), type = "closed")) {
  my_plot_theme() %+replace% 
    theme_noaxis(axis.line.x.bottom = element_line2(
    id = 1,
    xlength = xlength, arrow = arrow
  ), axis.line.y.left = element_line2(
    id = 2,
    ylength = ylength, arrow = arrow
  ), axis.title = element_text(hjust = 0.1,size = 7))
}

base_colors <- pal_npg()(10)
adaptive_palette <- function(n) {
  colorRampPalette(base_colors)(n)
}
