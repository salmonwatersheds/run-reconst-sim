library(here)

###############################################################################
# Load data
###############################################################################

CU_names <- read.csv("data/CentralCoastChum_CUs.csv")
CU_names <- CU_names[,'CU_name']

# Note: CU name is misspelled as Bella-Colla Dean in nccdbv2...I changed it
datAge <- read.csv("data/nccdbv2_age-composition-data.csv")

datAge <- subset(datAge, is.element(datAge$Sp.CU_Name, paste("CM", CU_names, sep="-")))

ages <- datAge[datAge$HasAgeData=="Yes", c('Age.2', "Age.3", "Age.4", "Age.5", "Age.6", "Age.7")]

# Check they sum to one
apply(ages, 1, sum, na.rm = TRUE)

# Calculate mean across all central coast chum CUs
# *This is the true age assumed in the model. *
apply(ages, 2, mean)
