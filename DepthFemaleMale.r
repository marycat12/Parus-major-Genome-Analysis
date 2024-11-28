rm(list=ls())
library (qqman)
library(ggplot2)
library(plyr)

c_m <- read.delim("MethylutionA-Male.merged-mapped.bam.amp-depth.txt", header = T)
c_m$chrom_name <- as.factor(c_m$X)
c_m$chrom_name <- revalue(c_m$chrom_name,
                          c(NC_031768.1="1",
                            NC_031769.1="2",
                            NC_031770.1="3",
                            NC_031771.1="4",
                            NC_031772.1="4A",
                            NC_031773.1="1A",
                            NC_031774.1="5",
                            NC_031775.1="6",
                            NC_031776.1="7",
                            NC_031777.1="8",
                            NC_031778.1="9",
                            NC_031779.1="10",
                            NC_031780.1="11",
                            NC_031781.1="12",
                            NC_031782.1="13",
                            NC_031783.1="14",
                            NC_031784.1="15",
                            NC_031785.1="17",
                            NC_031786.1="18",
                            NC_031787.1="19",
                            NC_031788.1="20",
                            NC_031789.1="21",
                            NC_031790.1="22",
                            NC_031791.1="23",
                            NC_031792.1="24",
                            NC_031793.1="25LG1",
                            NC_031794.1="25LG2",
                            NC_031795.1="26",
                            NC_031796.1="27",
                            NC_031797.1="28",
                            NC_031799.1="Z",
                            NC_031798.1="LGE22",
                            NC_040875.1="MT"
                          ))



c_f <- read.delim("MethylutionJ-Female.merged-mapped.bam.amp-depth.txt", header = T)
c_f$chrom_name <- as.factor(c_f$X)

c_f$chrom_name <- revalue(c_f$chrom_name,
                        c(NC_031768.1="1",
                        NC_031769.1="2",
                        NC_031770.1="3",
                        NC_031771.1="4",
                        NC_031772.1="4A",
                        NC_031773.1="1A",
                        NC_031774.1="5",
                        NC_031775.1="6",
                        NC_031776.1="7",
                        NC_031777.1="8",
                        NC_031778.1="9",
                        NC_031779.1="10",
                        NC_031780.1="11",
                        NC_031781.1="12",
                        NC_031782.1="13",
                        NC_031783.1="14",
                        NC_031784.1="15",
                        NC_031785.1="17",
                        NC_031786.1="18",
                        NC_031787.1="19",
                        NC_031788.1="20",
                        NC_031789.1="21",
                        NC_031790.1="22",
                        NC_031791.1="23",
                        NC_031792.1="24",
                        NC_031793.1="25LG1",
                        NC_031794.1="25LG2",
                        NC_031795.1="26",
                        NC_031796.1="27",
                        NC_031797.1="28",
                        NC_031799.1="Z",
                        NC_031798.1="LGE22",
                        NC_040875.1="MT"
                        ))

########### merge male and female dataframes #######

c_mf <- merge(c_m, c_f,
                 by=c("X",
                      "X.1",
                      "chrom_name"),
                 suffix = c("_m", "_f"))


######## look at all chromosomes and unplaced scaffolds ##########
# number of reads
corr_reads <- ggplot(c_mf,
                     aes(y = X.2_m,
                         x = X.2_f,
                         label=chrom_name)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point() +
    geom_text(hjust=0, vjust=0)
corr_reads



######## look only at all "major" chromosomes ##########

cc_mf <- c_mf[1:183146,]

cc_f <- c_f[1:996696,]

cc_m <- c_m[1:506226,]

# number of reads
corr_reads <- ggplot(cc_mf, aes(y = X.2_m, x = X.2_f, label = chrom_name)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  geom_text(hjust = -0.5, vjust = 1) +
  xlab("Amount of amplification for female") +  # Move xlab here
  ylab("Amount of amplification for male")  # Move ylab here

corr_reads


#Z Score
cc_mf$z.2_m=scale(cc_mf$X.2_m, center=TRUE, scale=TRUE)
cc_mf$z.2_f=scale(cc_mf$X.2_f, center=TRUE, scale=TRUE)

corr_reads <- ggplot(cc_mf,
                     aes(y = z.2_m,
                         x = z.2_f,
                         label=chrom_name)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  geom_text(hjust=-0.5, vjust=1) +
  xlab("Amount of amplification for female") +  # Move xlab here
  ylab("Amount of amplification for male")  # Move ylab here
corr_reads

############## sliding window analysis #########################

rm(list=ls())
library (qqman)
library(ggplot2)
library(plyr)
library(scales)

######## Male ##########

c_m <- read.delim("MethylutionA-Male.merged-mapped.bam.amp-depth.txt", header = T)

proper_chroms <- c("NC_031768.1",
                   "NC_031769.1",
                   "NC_031770.1",
                   "NC_031771.1",
                   "NC_031772.1",
                   "NC_031773.1",
                   "NC_031774.1",
                   "NC_031775.1",
                   "NC_031776.1",
                   "NC_031777.1",
                   "NC_031778.1",
                   "NC_031779.1",
                   "NC_031780.1",
                   "NC_031781.1",
                   "NC_031782.1",
                   "NC_031783.1",
                   "NC_031784.1",
                   "NC_031785.1",
                   "NC_031786.1",
                   "NC_031787.1",
                   "NC_031788.1",
                   "NC_031789.1",
                   "NC_031790.1",
                   "NC_031791.1",
                   "NC_031792.1",
                   "NC_031793.1",
                   "NC_031794.1",
                   "NC_031795.1",
                   "NC_031796.1",
                   "NC_031797.1",
                   "NC_031799.1",
                   "NC_031798.1",
                   "NC_040875.1")


c_m <- c_m[c_m$X %in% proper_chroms,]
c_m$chrom_num <- c_m$X



c_m$chrom_num  <- revalue(c_m$chrom_num,
                          c(NC_031768.1=1,
                            NC_031769.1=2,
                            NC_031770.1=3,
                            NC_031771.1=4,
                            NC_031772.1=4.5,
                            NC_031773.1=1.5,
                            NC_031774.1=5,
                            NC_031775.1=6,
                            NC_031776.1=7,
                            NC_031777.1=8,
                            NC_031778.1=9,
                            NC_031779.1=10,
                            NC_031780.1=11,
                            NC_031781.1=12,
                            NC_031782.1=13,
                            NC_031783.1=14,
                            NC_031784.1=15,
                            NC_031785.1=17,
                            NC_031786.1=18,
                            NC_031787.1=19,
                            NC_031788.1=20,
                            NC_031789.1=21,
                            NC_031790.1=22,
                            NC_031791.1=23,
                            NC_031792.1=24,
                            NC_031793.1=25.0,
                            NC_031794.1=25.5,
                            NC_031795.1=26,
                            NC_031796.1=27,
                            NC_031797.1=28,
                            NC_031799.1=29,
                            NC_031798.1=30,
                            NC_040875.1=31
                          ))
c_m$chrom_num <- as.numeric(c_m$chrom_num)
c_m$SNP <- NA

c_m$chromosome <- revalue(c_m$X,
                          c(NC_031768.1="1",
                            NC_031769.1="2",
                            NC_031770.1="3",
                            NC_031771.1="4",
                            NC_031772.1="4A",
                            NC_031773.1="1A",
                            NC_031774.1="5",
                            NC_031775.1="6",
                            NC_031776.1="7",
                            NC_031777.1="8",
                            NC_031778.1="9",
                            NC_031779.1="10",
                            NC_031780.1="11",
                            NC_031781.1="12",
                            NC_031782.1="13",
                            NC_031783.1="14",
                            NC_031784.1="15",
                            NC_031785.1="17",
                            NC_031786.1="18",
                            NC_031787.1="19",
                            NC_031788.1="20",
                            NC_031789.1="21",
                            NC_031790.1="22",
                            NC_031791.1="23",
                            NC_031792.1="24",
                            NC_031793.1="25LG1",
                            NC_031794.1="25LG2",
                            NC_031795.1="26",
                            NC_031796.1="27",
                            NC_031797.1="28",
                            NC_031799.1="Z",
                            NC_031798.1="LGE22",
                            NC_040875.1="MT"
                          ))

#Normal Axis
all_chroms_eggs_m_normal <- manhattan(c_m, chr = "chrom_num",
                                    bp = "X.1",
                                    p = "X.2",
                                    suggestiveline = 80,
                                    genomewideline = FALSE,
                                    cex = 0.7,
                                    logp = FALSE,
                                    col = c("blue", "red"),
                                    ylab = "amount of amplification",
                                    main = "Sample A (male)",
                                    chrlabs = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", 
                                                "8", "9", "10", "11", "12", "13", "14", "15", "17", 
                                                "18", "19", "20", "21", "22", "23", "24", "25LG1", 
                                                "26", "27", "28", "Z"),
                                    )

all_chroms_eggs_m_normal



#Logarithmische Achse
# Compute log10 transformation of the data
c_m$X.2_log <- log10(c_m$X.2)

# Define target regions (chromosome, start, end, and color)
target_regions <- data.frame(
  chrom_num = c(1, 3, 3, 28),
  start = c(24302500, 30299000, 30593000, 12500),
  end = c(24304500, 30301500, 30594500, 17500),
  color = c("green", "green", "green", "green")
)


# Custom Manhattan plot with target regions highlighted
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
                                         "18", "19", "20", "21", "22", "23", "24", "25LG1", "26", "27", "28", "Z"))


all_chroms_eggs   


###########TARGET REGIONS###########
###ALL
# Zielregionen definieren
target_regions <- data.frame(
  chrom_num = c(1, 3, 3, 28),
  start = c(24302500, 30299000, 30593000, 12500),
  end = c(24304500, 30301500, 30594500, 17500),
  color = c("green", "green", "green", "green")
)

# Farbcodierung der Daten basierend auf Zielregionen
c_m$color <- "blue"  # Standardfarbe für alle Punkte
for (i in 1:nrow(target_regions)) {
  region <- target_regions[i, ]
  c_m$color[c_m$chrom_num == region$chrom_num & c_m$X.1 >= region$start & c_m$X.1 <= region$end] <- region$color
}

# ggplot-Manhattan-Plot erstellen
ggplot(c_m, aes(x = chrom_num, y = X.2, color = color)) +
  geom_point(size = 1) +  # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) +  # Farben zuweisen
  scale_x_continuous(breaks = 1:31, labels = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24", "25LG1", "26", "27", "28", "Z", "MT")) +  # x-Achse korrekt beschriften
  labs(
    title = "Sample A (male) with Target Regions Highlighted",
    x = "Chromosome",
    y = "amount of amplification"
  ) +
  theme_minimal() +  # Minimales ggplot-Theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # X-Achse text drehen
    legend.position = "none"  # Legende ausblenden
  )

##ATTEMPT2
# Assign color based on SNP position and chromosome
c_m$color <- "blue"  # Default color is blue

# Apply green triangle shape for SNPs in the specified regions
c_m$shape <- 16  # Default shape is a circle (16) for all SNPs

# Apply triangle shape (17) for SNPs in the specified regions
c_m$shape[c_m$chrom_num == 3 & c_m$X.1 >= 30299000 & c_m$X.1 <= 30301500] <- 17
c_m$shape[c_m$chrom_num == 3 & c_m$X.1 >= 30593000 & c_m$X.1 <= 30594500] <- 17
c_m$shape[c_m$chrom_num == 1 & c_m$X.1 >= 24302500 & c_m$X.1 <= 24304500] <- 17
c_m$shape[c_m$chrom_num == 28 & c_m$X.1 >= 12500 & c_m$X.1 <= 17500] <- 17

# Assign green color to these triangles
c_m$color[c_m$shape == 17] <- "green"

# Now plot with the updated color column
all_chroms_eggs <- manhattan(c_m, chr = "chrom_num",
                             bp = "X.1",
                             p = "X.2_log",
                             suggestiveline = 80,
                             genomewideline = FALSE,
                             cex = 0.7,
                             logp = FALSE,
                             col = c_m$color,  # Use the updated color column
                             ylab = "log10 (amount of amplification)",
                             main = "Sample A (male)",
                             chrlabs = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", 
                                         "8", "9", "10", "11", "12", "13", "14", "15", "17", 
                                         "18", "19", "20", "21", "22", "23", "24", "25LG1", 
                                         "26", "27", "28", "Z"))

# Display the Manhattan plot with highlighted regions
all_chroms_eggs

table(c_m$color)

###Primerpair 1
# Daten für Chromosom 3 filtern
chrom1_data <- subset(c_m, chrom_num == 3)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 30299000 & chrom1_data$X.1 <= 30301500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample A (male) Chromosome 3 with target region of primer pair 1 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden

###Primerpair 2
# Daten für Chromosom 3 filtern
chrom1_data <- subset(c_m, chrom_num == 3)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 30593000 & chrom1_data$X.1 <= 30594500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample A (male) Chromosome 3 with target region of primer pair 2 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden

###Primerpair 6
# Daten für Chromosom 1 filtern
chrom1_data <- subset(c_m, chrom_num == 1)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 24302500 & chrom1_data$X.1 <= 24304500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample A (male) Chromosome 1 with target region of primer pair 6 in green",
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

######## Female ##########

c_f <- read.delim("MethylutionJ-Female.merged-mapped.bam.amp-depth.txt", header = T)

proper_chroms <- c("NC_031768.1",
                   "NC_031769.1",
                   "NC_031770.1",
                   "NC_031771.1",
                   "NC_031772.1",
                   "NC_031773.1",
                   "NC_031774.1",
                   "NC_031775.1",
                   "NC_031776.1",
                   "NC_031777.1",
                   "NC_031778.1",
                   "NC_031779.1",
                   "NC_031780.1",
                   "NC_031781.1",
                   "NC_031782.1",
                   "NC_031783.1",
                   "NC_031784.1",
                   "NC_031785.1",
                   "NC_031786.1",
                   "NC_031787.1",
                   "NC_031788.1",
                   "NC_031789.1",
                   "NC_031790.1",
                   "NC_031791.1",
                   "NC_031792.1",
                   "NC_031793.1",
                   "NC_031794.1",
                   "NC_031795.1",
                   "NC_031796.1",
                   "NC_031797.1",
                   "NC_031799.1",
                   "NC_031798.1",
                   "NC_040875.1")


c_f <- c_f[c_f$X %in% proper_chroms,]
c_f$chrom_num <- c_f$X



c_f$chrom_num  <- revalue(c_f$chrom_num,
                          c(NC_031768.1=1,
                            NC_031769.1=2,
                            NC_031770.1=3,
                            NC_031771.1=4,
                            NC_031772.1=4.5,
                            NC_031773.1=1.5,
                            NC_031774.1=5,
                            NC_031775.1=6,
                            NC_031776.1=7,
                            NC_031777.1=8,
                            NC_031778.1=9,
                            NC_031779.1=10,
                            NC_031780.1=11,
                            NC_031781.1=12,
                            NC_031782.1=13,
                            NC_031783.1=14,
                            NC_031784.1=15,
                            NC_031785.1=17,
                            NC_031786.1=18,
                            NC_031787.1=19,
                            NC_031788.1=20,
                            NC_031789.1=21,
                            NC_031790.1=22,
                            NC_031791.1=23,
                            NC_031792.1=24,
                            NC_031793.1=25.0,
                            NC_031794.1=25.5,
                            NC_031795.1=26,
                            NC_031796.1=27,
                            NC_031797.1=28,
                            NC_031799.1=29,
                            NC_031798.1=30,
                            NC_040875.1=31
                          ))
c_f$chrom_num <- as.numeric(c_f$chrom_num)
c_f$SNP <- NA

c_f$chromosome <- revalue(c_f$X,
                          c(NC_031768.1="1",
                            NC_031769.1="2",
                            NC_031770.1="3",
                            NC_031771.1="4",
                            NC_031772.1="4A",
                            NC_031773.1="1A",
                            NC_031774.1="5",
                            NC_031775.1="6",
                            NC_031776.1="7",
                            NC_031777.1="8",
                            NC_031778.1="9",
                            NC_031779.1="10",
                            NC_031780.1="11",
                            NC_031781.1="12",
                            NC_031782.1="13",
                            NC_031783.1="14",
                            NC_031784.1="15",
                            NC_031785.1="17",
                            NC_031786.1="18",
                            NC_031787.1="19",
                            NC_031788.1="20",
                            NC_031789.1="21",
                            NC_031790.1="22",
                            NC_031791.1="23",
                            NC_031792.1="24",
                            NC_031793.1="25LG1",
                            NC_031794.1="25LG2",
                            NC_031795.1="26",
                            NC_031796.1="27",
                            NC_031797.1="28",
                            NC_031799.1="Z",
                            NC_031798.1="LGE22",
                            NC_040875.1="MT"
                          ))

#Normal Axis
all_chroms_eggs_f_normal <- manhattan(c_f, chr = "chrom_num",
                                      bp = "X.1_K",
                                      p = "X.2",
                                      suggestiveline = 80,
                                      genomewideline = FALSE,
                                      cex = 0.7,
                                      logp = FALSE,
                                      col = c("blue", "red"),
                                      ylab = "amount of amplification",
                                      main = "Sample J (female)",
                                      chrlabs = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", 
                                                  "8", "9", "10", "11", "12", "13", "14", "15", "17", 
                                                  "18", "19", "20", "21", "22", "23", "24", "25LG1", 
                                                  "26", "27", "28", "Z", "LGE22"),
)

all_chroms_eggs_f_normal



#Logarithmische Achse
c_f$X.2_log <- log10(c_f$X.2)

# LOG
all_chroms_eggs_f_log <- manhattan(c_f, chr = "chrom_num",
                             bp = "X.1",
                             p = "X.2_log",
                             suggestiveline = 80,
                             genomewideline = FALSE,
                             cex = 0.7,
                             logp = FALSE,
                             col = c("blue", "red"),
                             ylab = "log10 (amount of amplification)",
                             main = "Sample J (female)",
                             chrlabs = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", 
                                         "8", "9", "10", "11", "12", "13", "14", "15", "17", 
                                         "18", "19", "20", "21", "22", "23", "24", "25LG1", 
                                         "26", "27", "28", "Z", "LGE22"))


all_chroms_eggs_f_log

###########TARGET REGIONS###########
###ALL
# Zielregionen definieren
target_regions <- data.frame(
  chrom_num = c(1, 3, 3, 28),
  start = c(24302500, 30299000, 30593000, 12500),
  end = c(24304500, 30301500, 30594500, 17500),
  color = c("green", "green", "green", "green")
)

# Farbcodierung der Daten basierend auf Zielregionen
c_f$color <- "blue"  # Standardfarbe für alle Punkte
for (i in 1:nrow(target_regions)) {
  region <- target_regions[i, ]
  c_f$color[c_f$chrom_num == region$chrom_num & c_f$X.1 >= region$start & c_f$X.1 <= region$end] <- region$color
}

# ggplot-Manhattan-Plot erstellen
ggplot(c_f, aes(x = chrom_num, y = X.2_log, color = color)) +
  geom_point(size = 1) +  # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) +  # Farben zuweisen
  scale_x_continuous(breaks = 1:31, labels = c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24", "25LG1", "26", "27", "28", "Z", "LGE22")) +  # x-Achse korrekt beschriften
  labs(
    title = "Sample J (female) with Target Regions Highlighted",
    x = "Chromosome",
    y = "amount of amplification"
  ) +
  theme_minimal() +  # Minimales ggplot-Theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # X-Achse text drehen
    legend.position = "none"  # Legende ausblenden
  )

###Primerpair 1
# Daten für Chromosom 3 filtern
chrom1_data <- subset(c_f, chrom_num == 3)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 30299000 & chrom1_data$X.1 <= 30301500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample J (female) Chromosome 3 with target region of primer pair 1 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden

###Primerpair 2
# Daten für Chromosom 3 filtern
chrom1_data <- subset(c_f, chrom_num == 3)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 30593000 & chrom1_data$X.1 <= 30594500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample J (female) Chromosome 3 with target region of primer pair 2 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden

###Primerpair 6
# Daten für Chromosom 1 filtern
chrom1_data <- subset(c_f, chrom_num == 1)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 24302500 & chrom1_data$X.1 <= 24304500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample J (female) Chromosome 1 with target region of primer pair 6 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden

### Primerpair 7
#  Daten für Chromosom 28 filtern
chrom1_data <- subset(c_f, chrom_num == 28)

# Daten mit einer Farbcodierung erweitern
chrom1_data$color <- ifelse(chrom1_data$X.1 >= 12500 & chrom1_data$X.1 <= 17500, "green", "blue")

# Manhattan-Plot mit ggplot2
ggplot(chrom1_data, aes(x = X.1, y = X.2_log, color = color)) +
  geom_point(size = 1) + # Punkte plotten
  scale_color_manual(values = c("blue" = "blue", "green" = "green")) + # Farben zuweisen
  labs(
    title = "Sample J (female) Chromosome 28 with target region of primer pair 7 in green",
    x = "Base Pair Position",
    y = "log10 (amount of amplification)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Legende ausblenden