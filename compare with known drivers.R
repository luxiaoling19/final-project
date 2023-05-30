library(readxl)
library(writexl)


COSMIC_known_drivers <- read_excel("C:/Users/92098/Desktop/driver genes data/COSMIC known drivers.xlsx")
colnames(COSMIC_known_drivers)[1] <- "gene_name"

##########################################
histological_driver <- read_excel("C:/Users/92098/Desktop/driver genes data/histological_driver.xlsx")

histological_driver_name <-na.omit(histological_driver[,1])
colnames(histological_driver_name)[1] <- "gene_name"
histological_driver_name <- histological_driver_name[!grepl("gene_name", histological_driver_name$gene_name), ]
histological_driver_name<-unique(histological_driver_name)
histological_known_driver <- merge(histological_driver_name, COSMIC_known_drivers, by="gene_name")

#save driver genes to excel document

folder_path <- "C:/Users/92098/Desktop/driver genes data"
excel_file <- file.path(folder_path, "histological_known_drivers.xlsx")
write_xlsx(histological_known_driver, excel_file)

#determine whether driver is known
result <- histological_driver_name$gene_name %in% COSMIC_known_drivers$gene_name
print(result)
table(result) #count number of TRUE and FALSE
#############################################
COSMIC_known_drivers_contain_lung <- COSMIC_known_drivers[grepl("gene_name", histological_driver_name$gene_name), ]

major_SBS <- read_excel("C:/Users/92098/Desktop/driver genes data/major SBS.xlsx")
major_SBS <- major_SBS[7]
colnames(major_SBS)[1] <- "gene_name"
major_SBS <- na.omit(major_SBS)
major_SBS_known_driver <- merge(major_SBS, COSMIC_known_drivers, by="gene_name")

folder_path <- "C:/Users/92098/Desktop/driver genes data"
excel_file <- file.path(folder_path, "major_SBS_known_drivers.xlsx")
write_xlsx(major_SBS_known_driver, excel_file)

#determine whether driver is known
result <- major_SBS$gene_name %in% COSMIC_known_drivers$gene_name
print(result)
table(result)
#############################################


SBS2_and_SBS13 <- read_excel("C:/Users/92098/Desktop/driver genes data/SBS2 & SBS13.xlsx", 
                         sheet = "SBS2和13同时>0.1")
SBS2_and_SBS13 <- SBS2_and_SBS13[,"gene_name"]
SBS2_and_SBS13_known_drivers <- merge(SBS2_and_SBS13, COSMIC_known_drivers, by="gene_name")

excel_file <- file.path(folder_path, "SBS2_and_SBS13_known_drivers.xlsx")
write_xlsx(SBS2_and_SBS13_known_drivers, excel_file)
#determine whether driver is known
result <- SBS2_and_SBS13$gene_name %in% COSMIC_known_drivers$gene_name
print(result)
table(result)
#############################################
SBS2_or_SBS13 <- read_excel("C:/Users/92098/Desktop/driver genes data/SBS2 & SBS13.xlsx", 
                         sheet = "SBS2或SBS13其一>0.1")
SBS2_or_SBS13 <- SBS2_or_SBS13[1:20,"gene_name"]
SBS2_or_SBS13_known_drivers <- merge(SBS2_or_SBS13, COSMIC_known_drivers, by="gene_name")

excel_file <- file.path(folder_path, "SBS2_or_SBS13_known_drivers.xlsx")
write_xlsx(SBS2_or_SBS13_known_drivers, excel_file)
#determine whether driver is known
result <- SBS2_or_SBS13$gene_name %in% COSMIC_known_drivers$gene_name
print(result)
table(result)
#############################################

SBS4_and_SBS92 <- read_excel("C:/Users/92098/Desktop/driver genes data/SBS4 & SBS 92.xlsx", 
                          sheet = "&")
SBS4_and_SBS92 <- na.omit(SBS4_and_SBS92[,"gene_name"])
SBS4_and_SBS92_known_drivers <- merge(SBS4_and_SBS92, COSMIC_known_drivers, by="gene_name")

excel_file <- file.path(folder_path, "SBS4_and_SBS92_known_drivers.xlsx")
write_xlsx(SBS4_and_SBS92_known_drivers, excel_file)
#determine whether driver is known
result <- SBS4_and_SBS92$gene_name %in% COSMIC_known_drivers$gene_name
print(result)
table(result)
#############################################
SBS4_or_SBS92 <- read_excel("C:/Users/92098/Desktop/driver genes data/SBS4 & SBS 92.xlsx", 
                             sheet = "或")
SBS4_or_SBS92 <- na.omit(SBS4_or_SBS92[1:30,"gene_name"])
SBS4_or_SBS92_known_drivers <- merge(SBS4_or_SBS92, COSMIC_known_drivers, by="gene_name")

excel_file <- file.path(folder_path, "SBS4_or_SBS92_known_drivers.xlsx")
write_xlsx(SBS4_or_SBS92_known_drivers, excel_file)
#determine whether driver is known
result <- SBS4_or_SBS92$gene_name %in% COSMIC_known_drivers$gene_name
print(result)
table(result)


columns_to_search <- c("Name", "Tumour Types(Somatic)", "Cancer Syndrome","Other Syndrome") 

COSMIC_known_drivers_contain_lung <- COSMIC_known_drivers[rowSums(sapply(COSMIC_known_drivers[columns_to_search], function(x) grepl("NSCLC", x))) > 0, ]
nrow(COSMIC_known_drivers_contain_lung)



