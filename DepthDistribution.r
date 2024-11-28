rm(list=ls())

###FEMALE
data <- read.table("MethylutionJ-Female.merged-mapped.bam.merge-amp-depth.bed", header = FALSE)
value_counts <- table(column_V4)

barplot(value_counts, 
        main = "Häufigkeitsdiagramm von V4", 
        xlab = "Werte in V4", 
        ylab = "Häufigkeit", 
        col = "blue", 
        border = "black", 
        ylim = c(0, 120)) # Y-Achse bis 120