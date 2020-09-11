rm(list = ls())
################
## 4CE NEURO ##
###############

# Set working directory where the files are
#setwd("./4CE/phase1.1/")

#############
# LIBRARIES #
#############
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)


#function to create a list with the files to analyze
fileList <- function(path, pattern, selection) {
  fileListInput <- list.files(path = path, pattern = pattern)
  #remove FICHOS and old files
  fileListInput  <-
    fileListInput[!grepl(paste(
      c("FICHOS", "VA.csv", "BCH.csv", "CHOP.csv", "RP401.csv") ,
      collapse = "|"
    ), x = fileListInput)]

  if (selection == "pediatric") {
    fileListInput  <-
      fileListInput[grepl(paste(c("PED", "UNCCH"), collapse = "|"), x = fileListInput)]
    fileListInput <-
      fileListInput[!grepl("APHPPED.csv", x = fileListInput)]
  } else if (selection == "adult") {
    fileListInput  <-
      fileListInput[!grepl(paste(c("PED"), collapse = "|"), x = fileListInput)]
  } else if (selection == "all") {
    fileListInput <-
      fileListInput[!grepl(paste(
        c(
          "APHPPED.csv",
          "APHPPEDHOSP.csv",
          "BCH.csv",
          "CHOP.csv",
          "FRBDXPED.csv",
          "MGBPED.csv",
          "NUHPED.csv",
          "NWUPED.csv",
          "RP401.csv",
          "UKFRPED.csv",
          "UMICHPED.csv",
          "UPITTPED.csv"
        ),
        collapse = "|"
      ), x = fileListInput)]
  }
  return(fileListInput)
}

# create a list with all the files
# selection argument options: all, pediatric, adult
fileListDiag <-
  fileList(path = "./phase1.1/latest/",
           pattern = "Diag",
           selection = "all")
fileListDemog <-
  fileList(path = "./phase1.1/latest/",
           pattern = "Demog",
           selection = "all")

#put together all diagnosis data
icd9 <- 0
icd10 <- 0
icd9and10 <- 0

for (i in 1:length(fileListDiag)) {
  print(fileListDiag[i])
  selection <-
    read.delim(
      paste0("./phase1.1/latest/", fileListDiag[i]),
      sep = ",",
      colClasses = "character"
    )
  colnames(selection) <- tolower(colnames(selection))
  print(paste0(unique(selection$icd_version), collapse = ","))

  if (length(unique(selection$icd_version)) == 1) {
    if (unique(selection$icd_version) == "9") {
      icd9 <- icd9 + 1
    } else if (unique(selection$icd_version) == "10") {
      icd10 <- icd10 + 1
    }
  } else if (length(unique(selection$icd_version)) == 2) {
    icd9and10 <- icd9and10 + 1
  }

  if (i == 1) {
    allDiagnosis <- selection
  } else{
    allDiagnosis <- rbind(allDiagnosis, selection)
  }
}

#put together all demographic data
for (i in 1:length(fileListDemog)) {
  print(i)
  selection <-
    read.delim(
      paste0("./phase1.1/latest/", fileListDemog[i]),
      sep = ",",
      colClasses = "character"
    )
  colnames(selection) <- tolower(colnames(selection))
  selection <- selection[selection$sex == "all" &
                           selection$age_group == "all" &
                           selection$race == "all",]

  if (i == 1) {
    allDemographics <- selection
  } else {
    allDemographics <- rbind(allDemographics, selection)
  }
}


#estimate the total never severe per site and change the column names to know they are refered to total counts per site
allDemographics[allDemographics < 0] <- NA
allDemographics$num_patients_never_severe <-
  as.numeric(allDemographics$num_patients_all) - as.numeric(allDemographics$num_patients_ever_severe)
allDemographics <-
  allDemographics[, c(
    "siteid",
    "num_patients_all",
    "num_patients_ever_severe",
    "num_patients_never_severe"
  )]
colnames(allDemographics) <-
  c("siteid",
    "totalPatients",
    "totalEverSever",
    "totalNeverSevere")

###put all the information together merging all dataframes by site id
finalDataSet <- merge(allDiagnosis, allDemographics, by = "siteid")
finalDataSet$siteid <- toupper(finalDataSet$siteid)


########################################
## ICSM: merge together all ICSM data ##
########################################
ICSMdata <-
  finalDataSet[finalDataSet$siteid %in% c("ICSM1", "ICSM20", "ICSM5"),-1]
#transform -99 and -999 in NA
ICSMdata[, c(3:9)] <- sapply(ICSMdata[, c(3:9)], as.numeric)
ICSMdata[ICSMdata < 0] <- NA

ICSM <- ICSMdata %>%
  group_by(icd_code_3chars, icd_version) %>%
  #summarise_all( sum ) %>%
  dplyr::summarise_all( ~ {
    sum(.x, na.rm = any(!is.na(.x)))
  }) %>%
  mutate(siteid = "ICSM")

#remove the individual ICSM sets and add the aggregated one
ICSM <- as.data.frame(ICSM)
finalDataSet <-
  finalDataSet[!finalDataSet$siteid %in% c("ICSM1", "ICSM20", "ICSM5"),]
finalDataSet <- rbind(finalDataSet, ICSM)

#how many sites?
length(unique(finalDataSet$siteid))

###read the neuro-related ICD codes
neuroICD10 <-
  read_excel("./neuro/COVID_neuro_ICD_20200826.xlsx", sheet = 1)
neuroICD10 <- as.data.frame(neuroICD10)
neuroICD10 <- neuroICD10[neuroICD10[, 2] == 1, c(1, 3, 4)]
neuroICD10$ICD_version <- 10
colnames(neuroICD10)[2] <- "ICD_code"


#neuroICD9  <- read.delim( "./neuro/neurology10to9.txt")
#neuroICD9$digits <- apply(neuroICD9[1], 1, function(x) nchar(x))
#neuroICD9 <- unique( neuroICD9[ neuroICD9$digits == 3, c("ICD9_CODE", "ICD10_CODE" ) ] )

#neuroICD9 <- merge( neuroICD9, neuroICD10, by.x = "ICD10_CODE", by.y = "ICD_code")
#neuroICD9$ICD_version <- 9
#neuroICD9 <- neuroICD9[ , c(3,2,4,5) ]
#colnames(neuroICD9 )[2] <- "ICD_code"

#all neuro ICD codes
#allICDneuroCodes <- rbind( neuroICD10, neuroICD9 )
allICDneuroCodes <- neuroICD10

###select from the 4CE data the subset of codes
selection <-
  finalDataSet[finalDataSet$icd_code_3chars %in% allICDneuroCodes$ICD_code,]
selection <-
  merge(selection, allICDneuroCodes, by.x = "icd_code_3chars", by.y = "ICD_code")

##de-identify
mapping <-
  read.delim("./phase1.1/mappingFiles/SiteID_Map.csv", sep = ",")
mapping$Acronym <- toupper(mapping$Acronym)

###merge
selection <-
  merge(selection, mapping, by.x = "siteid", by.y = "Acronym")


###organize the columns
selection <-
  selection[, c(
    "Neurological Disease Category",
    "Description",
    "icd_code_3chars",
    "Anonymous.Site.ID",
    "Country",
    "totalPatients",
    "totalEverSever",
    "num_patients_all_before_admission",
    "num_patients_all_since_admission",
    "num_patients_ever_severe_before_admission",
    "num_patients_ever_severe_since_admission"
  )]
selection <-
  selection[order(selection$icd_code_3chars, decreasing = TRUE),]


##################################
## Obfuscated and unknow values ##
##################################
#transform -99 and -999 in NA
selection[, c(1:5)] <- sapply(selection[, c(1:5)], as.character)
selection[, c(6:11)] <- sapply(selection[, c(6:11)], as.numeric)
selection[selection < 0] <- NA

##### prevalence for each code for each site (before and after admission)
prevalenceBySite <- selection[, 1:9]
prevalenceBySite$BeforeAdmission <-
  round(
    100 * (
      prevalenceBySite$num_patients_all_before_admission / prevalenceBySite$totalPatients
    ),
    2
  )
prevalenceBySite$AfterAdmission <-
  round(
    100 * (
      prevalenceBySite$num_patients_all_since_admission / prevalenceBySite$totalPatients
    ),
    2
  )

#### Heatmap plot
head(prevalenceBySite)
prevalenceBySite$description <-
  paste0(prevalenceBySite$Description,
         "-",
         prevalenceBySite$icd_code_3chars)
toplotBA <-
  prevalenceBySite[, c(
    "Anonymous.Site.ID",
    "icd_code_3chars" ,
    "num_patients_all_before_admission",
    "BeforeAdmission",
    "Neurological Disease Category"
  )]
colnames(toplotBA) <-
  c("Site",
    "Description",
    "PatientCount",
    "Percentage",
    "DiseaseCategory")
toplotBA$When <- "Before Admission"
toplotAA <-
  prevalenceBySite[, c(
    "Anonymous.Site.ID",
    "icd_code_3chars" ,
    "num_patients_all_since_admission",
    "AfterAdmission",
    "Neurological Disease Category"
  )]
colnames(toplotAA) <-
  c("Site",
    "Description",
    "PatientCount",
    "Percentage",
    "DiseaseCategory")
toplotAA$When <- "After Admission"

toplot <- rbind(toplotAA, toplotBA)

#before admission
ggplot(data =  toplotBA, aes(x = Description, y = Site, fill = Percentage)) +
  geom_tile(aes(fill = Percentage), colour = "white") +
  geom_text(aes(label = PatientCount), col = 'black', cex = 2) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 6
    ),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    text = ggplot2::element_text(size = 8),
    axis.title = ggplot2::element_text(size = 6)
  ) +
  scale_fill_gradient(name = "Percentage",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_grid(~ DiseaseCategory, scales = "free")


#after admission
ggplot(data =  toplotAA, aes(x = Description, y = Site, fill = Percentage)) +
  geom_tile(aes(fill = Percentage), colour = "white") +
  geom_text(aes(label = PatientCount), col = 'black', cex = 2) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 6
    ),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    text = ggplot2::element_text(size = 8),
    axis.title = ggplot2::element_text(size = 6)
  ) +
  scale_fill_gradient(name = "Percentage",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_grid(~ DiseaseCategory, scales = "free")


##### prevalence for each code by country (before and after admission)
byCountry <- selection[, c(1:3, 5:9)]
prevalenceByCountry <- byCountry %>%
  group_by(Country,
           `Neurological Disease Category`,
           Description,
           icd_code_3chars) %>%
  dplyr::summarise_all( ~ {
    sum(.x, na.rm = any(!is.na(.x)))
  })

prevalenceByCountry <- as.data.frame(prevalenceByCountry)
prevalenceByCountry$BeforeAdmission <-
  round(
    100 * (
      prevalenceByCountry$num_patients_all_before_admission / prevalenceByCountry$totalPatients
    ),
    2
  )
prevalenceByCountry$AfterAdmission <-
  round(
    100 * (
      prevalenceByCountry$num_patients_all_since_admission / prevalenceByCountry$totalPatients
    ),
    2
  )

#### Heatmap plot

head(prevalenceByCountry)
prevalenceByCountry$description <-
  paste0(prevalenceByCountry$Description,
         "-",
         prevalenceByCountry$icd_code_3chars)
toplotBA <-
  prevalenceByCountry[, c(
    "icd_code_3chars" ,
    "Country",
    "num_patients_all_before_admission",
    "BeforeAdmission"
  )]
colnames(toplotBA) <-
  c("Description", "Country", "PatientCount", "Percentage")
toplotBA$When <- "Before Admission"
toplotAA <-
  prevalenceByCountry[, c(
    "icd_code_3chars" ,
    "Country",
    "num_patients_all_since_admission",
    "AfterAdmission"
  )]
colnames(toplotAA) <-
  c("Description", "Country", "PatientCount", "Percentage")
toplotAA$When <- "After Admission"

toplot <- rbind(toplotAA, toplotBA)

#heatmap
ggplot(data =  toplot, aes(x = When, y = Description, fill = Percentage)) +
  geom_tile(aes(fill = Percentage), colour = "white") +
  geom_text(aes(label = PatientCount), col = 'black', cex = 3) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_fill_gradient(name = "Percentage",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_grid(~ Country, scales = "free")

#dot plot
ggplot(data =  toplot,
       aes(
         x = Description,
         y = Percentage,
         color = When,
         group = Description
       )) +
  geom_point(alpha = 0.5) +
  geom_line(color = "black") +
  coord_flip() +
  facet_grid(~ Country) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )


##### prevalence for each code all (before and after admission)
byAll <- selection[, c(1:3, 6:9)]

prevalenceByAll <- byAll %>%
  group_by(`Neurological Disease Category`,
           Description, icd_code_3chars) %>%
  dplyr::summarise_all( ~ {
    sum(.x, na.rm = any(!is.na(.x)))
  })

prevalenceByAll <- as.data.frame(prevalenceByAll)

prevalenceByAll$BeforeAdmission <-
  round(
    100 * (
      prevalenceByAll$num_patients_all_before_admission / prevalenceByAll$totalPatients
    ),
    2
  )
prevalenceByAll$AfterAdmission <-
  round(
    100 * (
      prevalenceByAll$num_patients_all_since_admission / prevalenceByAll$totalPatients
    ),
    2
  )


#### Heatmap plot
head(prevalenceByAll)
prevalenceByAll$description <-
  paste0(prevalenceByAll$Description,
         "-",
         prevalenceByAll$icd_code_3chars)
toplotBA <-
  prevalenceByAll[, c("description" ,
                      "num_patients_all_before_admission",
                      "BeforeAdmission")]
colnames(toplotBA) <- c("Description", "PatientCount", "Percentage")
toplotBA$When <- "Before Admission"
toplotAA <-
  prevalenceByAll[, c("description" ,
                      "num_patients_all_since_admission",
                      "AfterAdmission")]
colnames(toplotAA) <- c("Description", "PatientCount", "Percentage")
toplotAA$When <- "After Admission"

toplot <- rbind(toplotAA, toplotBA)

##heatmap
ggplot(data =  toplot, aes(x = When, y = Description, fill = Percentage)) +
  geom_tile(aes(fill = Percentage), colour = "white") +
  geom_text(aes(label = PatientCount), col = 'black', cex = 3) +
  scale_fill_gradient(name = "Percentage",
                      low = "#FFFFFF",
                      high = "#012345")

##dot plot
ggplot(data =  toplot,
       aes(
         x = Description,
         y = Percentage,
         color = When,
         group = Description
       )) +
  geom_point() +
  geom_line(color = "black") +
  coord_flip()
