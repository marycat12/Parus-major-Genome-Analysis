rm(list=ls())  # Clear workspace

### FEMALE
# Load data
data <- read.table("MethylutionJ-Female.merged-mapped.bam.merge-amp-depth.bed", header = FALSE)

# Calculate the absolute length difference between V2 and V3
data$length_difference <- abs(data$V3 - data$V2)

# Create a frequency table of the length differences
value_counts <- table(data$length_difference)

# Create the bar plot with logarithmic y-axis
bar_positions <- barplot(value_counts, 
                         main = "Length Difference Distribution in Female Sample (J)", 
                         xlab = "Length Difference", 
                         ylab = "Log (Frequency)", 
                         col = "blue", 
                         border = "black", 
                         log = "y")  # Logarithmic scale

### MALE
# Load data
data <- read.table("MethylutionA-Male.merged-mapped.bam.merge-amp-depth.bed", header = FALSE)

# Calculate the absolute length difference between V2 and V3
data$length_difference <- abs(data$V3 - data$V2)

# Create a frequency table of the length differences
value_counts <- table(data$length_difference)

# Create the bar plot with logarithmic y-axis
bar_positions <- barplot(value_counts, 
                         main = "Length Difference Distribution in Male Sample (A)", 
                         xlab = "Length Difference", 
                         ylab = "Log (Frequency)", 
                         col = "blue", 
                         border = "black", 
                         log = "y")  # Logarithmic scale
