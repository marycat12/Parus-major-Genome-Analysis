rm(list=ls())

###FEMALE
data <- read.table("MethylutionJ-Female.merged-mapped.bam.merge-amp-depth.bed", header = FALSE)
value_counts <- table(data$V4)

#log
bar_positions <- barplot(value_counts, 
                         main = "Depth Distribution in Female Sample (J)", 
                         xlab = "Depth", 
                         ylab = "Log (Frequency)", 
                         col = "blue", 
                         border = "black", 
                         log = "y") # Logarithmische Skala

# Group depth values into 1-unit bins (adjusting the cut function)
data$V4_binned <- cut(data$V4, 
                      breaks = seq(floor(min(data$V4)), ceiling(max(data$V4)), by = 1), 
                      include.lowest = TRUE, right = TRUE)  # right = TRUE for including upper boundary

# Count the occurrences in each bin
value_counts <- table(data$V4_binned)

# Remove bins with zero counts
value_counts <- value_counts[value_counts > 0]

# Create the bar plot with logarithmic y-axis
bar_positions <- barplot(value_counts, 
                         main = "Depth Distribution in Female Sample (J)", 
                         xlab = "Depth (binned by 1)", 
                         ylab = "Log (Frequency)", 
                         col = "blue", 
                         border = "black", 
                         log = "y", 
                         las = 2, 
                         names.arg = sub("^\\((\\d+),.*", "\\1", names(value_counts)))  # Extract the lower boundary of each bin

#normal
barplot(value_counts, 
        main = "Depth Distribution in Female Sample (J)", 
        xlab = "Depth", 
        ylab = "Frequency", 
        col = "blue", 
        border = "black", 
        ylim = c(0, 120)) # Y-Achse bis 120

###MALE
data <- read.table("MethylutionA-Male.merged-mapped.bam.merge-amp-depth.bed", header = FALSE)
value_counts <- table(data$V4)

#log
bar_positions <- barplot(value_counts, 
                         main = "Depth Distribution in Male Sample (A)", 
                         xlab = "Depth", 
                         ylab = "Log (Frequency)", 
                         col = "blue", 
                         border = "black", 
                         log = "y") # Logarithmische Skala

# Group depth values into 1-unit bins (adjusting the cut function)
data$V4_binned <- cut(data$V4, 
                      breaks = seq(floor(min(data$V4)), ceiling(max(data$V4)), by = 1), 
                      include.lowest = TRUE, right = TRUE)  # right = TRUE for including upper boundary

# Count the occurrences in each bin
value_counts <- table(data$V4_binned)

# Remove bins with zero counts
value_counts <- value_counts[value_counts > 0]

# Create the bar plot with logarithmic y-axis
bar_positions <- barplot(value_counts, 
                         main = "Depth Distribution in Male Sample (A)", 
                         xlab = "Depth (binned by 1)", 
                         ylab = "Log (Frequency)", 
                         col = "blue", 
                         border = "black", 
                         log = "y", 
                         las = 2, 
                         names.arg = sub("^\\((\\d+),.*", "\\1", names(value_counts)))  # Extract the lower boundary of each bin

#normal
barplot(value_counts, 
        main = "Depth Distribution in Female Sample (J)", 
        xlab = "Depth", 
        ylab = "Frequency", 
        col = "blue", 
        border = "black", 
        ylim = c(0, 1200)) # Y-Achse bis 120
