library(here)

###############################################################################
# Load data
###############################################################################

CU_names <- read.csv("data/CentralCoastChum_CUs.csv")
CU_names <- CU_names[,1]

# Note: CU name is misspelled as Bella-Colla Dean in nccdbv2...I changed it
datAge <- read.csv("data/nccdbv2_age-composition-data.csv")

datAge <- subset(datAge, is.element(datAge$Sp.CU_Name, paste("CM", CU_names, sep="-")))

ages <- datAge[datAge$HasAgeData=="Yes", c('Age.2', "Age.3", "Age.4", "Age.5", "Age.6", "Age.7")]
# Check they sum to one
apply(ages, 1, sum, na.rm = TRUE)

# write.csv(cbind(c(as.character(datAge$Sp.CU_Name[datAge$HasAgeData=="Yes"]), "Average"), rbind(ages,apply(ages, 2, mean))), file="ages.csv")
sum(apply(ages, 2, mean), na.rm=TRUE)

sum(round(apply(ages, 2, mean), 2), na.rm=TRUE)
