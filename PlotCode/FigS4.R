tht_mat = readRDS("Output/THT/Tht_mat_1.RDS")
trackers = readRDS("Output/POST_TRACK/Post_track_1.RDS")

peaks = get.only.max.peaks(trackers$Tot_tht_sig, midline = 0.3)

tick = peaks[length(peaks) - 1]

tht_vec = tht_mat[tick, which(tht_mat[tick, ] != -1)]

rescale = function(x, vector = tht_vec) {
  (x - min(vector))/(max(vector) - min(vector))
}

cutoff = median(tht_vec) + sig_cutoff
cutoff = rescale(cutoff, tht_vec)

tht_vec = rescale(tht_vec) %>% as.data.frame()
colnames(tht_vec) = "ThT"

dat = with(density(tht_vec$ThT, from = 0,
                   to = 1), data.frame(x, y))

p_mem = ggplot(data = dat, mapping = aes(x = x, y = y)) +
  geom_area(mapping = aes(x = ifelse(x <= cutoff + 0.0004, x, 0),
                          fill = "Nonsignaling"), alpha = 0.9) +
  geom_area(mapping = aes(x = ifelse(x >= cutoff, x, -10), fill = "Signaling"),
            alpha = 0.9) +
  geom_vline(aes(xintercept = cutoff),
             size = 1, linetype = "dotted") +
  geom_line(color = "black") +
  scale_fill_manual(values = c("Signaling" = "slategray",
                                "Nonsignaling" = "lightgray")) +
  xlim(min(tht_vec), max(tht_vec)) +
  labs(title = "",
       x = "Scaled ThT",
       y = "Density",
       fill = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.position = c(0.86, 0.877)
  )

p_mem

ggsave(paste0(return_folder, "FigS4.tiff"), p_mem,
       width = 5.2, height = 4, units = "in",
       dpi = 300)
