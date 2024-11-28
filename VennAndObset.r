rm(list=ls())
library(ggVennDiagram)

MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$Region=paste(MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V1, MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V2, MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V3, sep=":")
MethylutionA.Male.merged.mapped.bam.merge.amp.depth$Region=paste(MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V1, MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V2, MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V3, sep=":")

#VENN DIAGRAMM
#ggVennDiagram(list(Male=MethylutionA.Male.merged.mapped.bam.merge.amp.depth$Region, Female=MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$Region))

#UPSET PLOT
ggVennDiagram(list(Male=MethylutionA.Male.merged.mapped.bam.merge.amp.depth$Region, Female=MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$Region), force_upset=TRUE)

#FEMALE Round on 100
# Round V2 and V3 to the nearest 100
MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V2 <- round(MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V2 / 100) * 100
MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V3 <- round(MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V3 / 100) * 100

# Recreate the Region column with the rounded V2 and V3
MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$Region <- paste(
  MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V1, 
  MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V2, 
  MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$V3, 
  sep = ":"
)

#MALE round on 100
# For Male dataset
MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V2 <- round(MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V2 / 100) * 100
MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V3 <- round(MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V3 / 100) * 100

# Recreate the Region column for Male
MethylutionA.Male.merged.mapped.bam.merge.amp.depth$Region <- paste(
  MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V1, 
  MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V2, 
  MethylutionA.Male.merged.mapped.bam.merge.amp.depth$V3, 
  sep = ":"
)

#VENN DIAGRAMM
#ggVennDiagram(list(Male=MethylutionA.Male.merged.mapped.bam.merge.amp.depth$Region, Female=MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$Region))

#UPSET PLOT
ggVennDiagram(list(Male=MethylutionA.Male.merged.mapped.bam.merge.amp.depth$Region, Female=MethylutionJ.Female.merged.mapped.bam.merge.amp.depth$Region), force_upset=TRUE)