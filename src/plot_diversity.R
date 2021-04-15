#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
    library(wesanderson)
    library(Rmisc)
    library(docopt)
})

#require(ggplot2); require(cowplot); 
#require(wesanderson); require(Rmisc); require(docopt)

'Plot Inverse Simpson diverity from VDJTOOLS CalcDiversityStats analysis
Usage:
   plot_diversity.R [-i <input> -c <chain> -t <threshold> -o <output>]

Arguments:
   -i Resampled diversity analysis input (tsv)
   -c Ig type [IGH, IGK, IGL]
   -t Read cutoff threshold
   -o Output plot

' -> doc

opts <- docopt(doc)

#args <- c("analysis/IgH/PlotFancyVJUsage/D02220_T2_IgH.fancyvj.wt.txt", "analysis/IgH/PlotFancyVJUsage/D02220_T2_IgH.fancyvj.wt.pdf")
#args<-commandArgs(TRUE)

file_in  <- opts$i
file_out <- opts$o
chain <- opts$c
cutoff <- opts$t

resamp <- read.csv(file_in, sep = "\t")

# swap D02220 time points
resamp[resamp$sample_id==paste("D02220_T1_",chain,sep=""),"time"] <- "T2"
resamp[resamp$sample_id==paste("D02220_T2_",chain,sep=""),"time"] <- "T1"

# print message to screen
cat("Plotting...", file_in, "\n", sep="")

# subset data
resamp_cases <- resamp[resamp$case_status == "1",]
resamp_cons <- resamp[resamp$case_status == "0",]

# calculate inverse simpson index CI
upper <- 0
lower <- 0
if (dim(resamp_cons)[1] != 0) {
    cis <- CI(resamp_cons$inverseSimpsonIndex_mean,ci=0.95)
    upper <- cis[1]
    lower <- cis[3]
}

# set plot color - PROBABLY CLEANER WAY TO DO THIS...
if (chain == "IGH") {
    color = "Darjeeling2"
} else if (chain == "IGK") {
    color = "GrandBudapest2"
} else {
    color = "Cavalcanti1"
}

# case plot
cases <- ggplot(resamp_cases, 
                aes(x = factor(dog_id), y = inverseSimpsonIndex_mean, fill = factor(time))) +
  scale_fill_manual(values = wes_palette(n = 4, name = color)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_text(size = 3, position = position_dodge(width=1), aes(y=inverseSimpsonIndex_mean+5, label=diversity, hjust=0), angle=45) +
  scale_x_discrete(limits = unique(resamp[resamp$case_status == "1" & 
                              resamp$diagnosis == "HGL",]$dog_id)) +
  labs(NULL) +
  ylim(0, max(resamp$inverseSimpsonIndex_mean) + 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# control plot
cons <- ggplot(resamp_cons, 
               aes(x = factor(dog_id), y = inverseSimpsonIndex_mean, fill = factor(time))) +
  scale_fill_manual(values = wes_palette(n = 4, name = color)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.5)) +
  geom_text(size = 3, position = position_dodge(width=1), aes(y=inverseSimpsonIndex_mean+5, label=diversity, hjust=0), angle=45) +
  scale_x_discrete(limits = unique(resamp[resamp$case_status == "0",]$dog_id)) +
  annotate("rect", 
           xmin = -Inf, xmax = Inf, 
           ymin = max(resamp$inverseSimpsonIndex_mean) + 5, ymax = upper, 
           fill = "grey", alpha = .5, color = NA) +
  geom_hline(yintercept = lower, linetype="dashed") +
  ylab(NULL) +
  ylim(0, max(resamp$inverseSimpsonIndex_mean) + 5) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

################################################################################
# plotting is dependent on if the cases and controls contain data at this 
# threshold
################################################################################
if ((dim(resamp_cons)[1] != 0) && (dim(resamp_cases) != 0)) {
    # combine case and control plots
    prow <- plot_grid(
      cases + theme(legend.position = "none"),
      cons + theme(legend.position = "none"),
      align = "vh",
      labels = c("Cases","Controls"),
      hjust = 0,
      nrow = 1
    )

    # generate plot, legend, and title
    legend <- get_legend(cases + theme(legend.box.margin = margin(0, 0, 0, 12)))
    p_reads <- plot_grid(prow, legend, rel_widths = c(3, .4))

    title <- ggdraw() + 
      draw_label(
        paste(chain," - cutoff ", cutoff, " (downsampled to ",resamp$resample_reads[1],")",sep = ""),
        fontface = 'bold',
        size = 15,
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 150)
      )

    plot_grid(
      title, p_reads,
      ncol = 1,
      # rel_heights values control vertical title margins
      rel_heights = c(0.1, 1)
    )

    ggsave(file_out, width = 12, height = 7)
    cat("Done. Output located at ", file_out, "\n")
# if there is case data, print cases only

} else if (dim(resamp_cases)[1] != 0) {
    # plot cases only
    prow <- plot_grid(
    cases + theme(legend.position = "none"),
    align = "vh",
    labels = c("Cases","Controls"),
    hjust = 0,
    nrow = 1
    )

    # generate plot, legend, and title
    legend <- get_legend(cases + theme(legend.box.margin = margin(0, 0, 0, 12)))
    p_reads <- plot_grid(prow, legend, rel_widths = c(3, .4))

    title <- ggdraw() + 
    draw_label(
      paste(chain," - cutoff ", cutoff, " (downsampled to ",resamp$resample_reads[1],")",sep = ""),
      fontface = 'bold',
      size = 15,
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 150)
    )

    plot_grid(
    title, p_reads,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
    )
    ggsave(file_out, width = 12, height = 7)
    cat("WARNING: only case data included at", cutoff, "\n")
    cat("Output located at ", file_out, "\n")

} else if (dim(resamp_cons)[1] != 0) {
    # plot controls only
    prow <- plot_grid(
    cons + theme(legend.position = "none"),
    align = "vh",
    labels = c("Cases","Controls"),
    hjust = 0,
    nrow = 1
    )

    # generate plot, legend, and title
    legend <- get_legend(cons + theme(legend.box.margin = margin(0, 0, 0, 12)))
    p_reads <- plot_grid(prow, legend, rel_widths = c(3, .4))

    title <- ggdraw() + 
    draw_label(
      paste(chain," - cutoff ", cutoff, " (downsampled to ",resamp$resample_reads[1],")",sep = ""),
      fontface = 'bold',
      size = 15,
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 150)
    )

    plot_grid(
    title, p_reads,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
    )
    ggsave(file_out, width = 12, height = 7)
    cat("WARNING: only control data included at", cutoff, "\n")
    cat("Output located at ", file_out, "\n")

} else {
    file.create(file_out)
    cat("WARNING: no case or control data included at", cutoff, "\n")
    cat("WARNING: empty output located at ", file_out, "\n")
}
