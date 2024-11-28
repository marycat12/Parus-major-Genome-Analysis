# Target regions
target_regions <- data.frame(
  chrom_num = c(1, 3, 3, 28),
  start = c(24302500, 30299000, 30593000, 12500),
  end = c(24304500, 30301500, 30594500, 17500),
  color = c("green", "green", "green", "green")
)

all_chroms_eggs <- manhattan(c_m, chr = "chrom_num",
                             bp = "X.1",
                             p = "X.2_log",
                             suggestiveline = 80,
                             genomewideline = FALSE,
                             cex = 0.7,
                             logp = FALSE,
                             col = c("blue", "red"),
                             ylab = "log10 (amount of amplification)",
                             main = "Sample A (male)",
                             chrlabs = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", 
                                         "8", "9", "10", "11", "12", "13", "14", "15", "17", 
                                         "18", "19", "20", "21", "22", "23", "24", "25LG1", 
                                         "26", "27", "28", "Z"))

# Display the initial Manhattan plot
all_chroms_eggs


####################################
###ALL
c_m$color1 <- ifelse(c_m$chrom_num == 3 & c_m$X.1 >= 30299000 & c_m$X.1 <= 30301500, "green", "blue")
c_m$color2 <- ifelse(c_m$chrom_num == 3 & c_m$X.1 >= 30593000 & c_m$X.1 <= 30594500, "green", "blue")
c_m$color6 <- ifelse(c_m$chrom_num == 1 & c_m$X.1 >= 24302500 & c_m$X.1 <= 24304500, "green", "blue")
c_m$color7 <- ifelse(c_m$chrom_num == 28 & c_m$X.1 >= 12500 & c_m$X.1 <= 17500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(c_m, aes(x = X.1, y = X.2_log, color1 = color1, color2 = color2, color6 = color6, color7 = color7)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample A (male) with targeted regions in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden


###Primerpair 7
# Daten für Chromosom 28 filtern
chrom1_data <- subset(c_m, chrom_num == 28)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 12500 & chrom1_data$X.1 <= 17500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample A (male) Chromosome 28 with target region of primer pair 7 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden
