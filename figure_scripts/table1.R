#############################################################
# Generate Table1 summary statistics for MMTM patients.  
# March 28, 2024
#############################################################

library(table1)

# Read the CSV file
data <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/mmtm_clinical_metadata.csv")


labels <- c(
  age = "Age",
                 bmi1 = "BMI",
                 ga_exact = "Exact Gestational Age",
                 nullip = "Nulliparity",
                 prior_ptb = "Prior Preterm Birth",
                 prior_sptb = "Prior Spontaneous Preterm Birth",
                 prior_term = "Prior Term Birth",
                 v1_ga = "GA at Visit 1",
                 v2_ga = "GA at Visit 2")

label(data$age)<- "Age"
label(data$bmi1)<- "BMI"
label(data$ga_exact)<- "Gestational Age at Delivery"
label(data$nullip)<- "Nulliparity"
label(data$prior_ptb)<- "Prior Preterm Births"
label(data$prior_sptb)<- "Prior Spontaneous Preterm Births"
label(data$prior_term)<- "Prior Term Births"
label(data$v1_ga)<- "Visit 1 Gestational Age"
label(data$v2_ga)<- "Vist 2 Gestational Age"

# Generate summary statistics using table1 function
summary_stats <- table1(~age + bmi1 + ga_exact + nullip + prior_ptb + prior_sptb + prior_term + v1_ga + v2_ga, data = data)

# Print the summary statistics
print(summary_stats)


# Generate summary statistics using table1 function
summary_stats <- table1(~nullip + prior_ptb + prior_sptb + prior_term, data = data, cols=4)

# Print the summary statistics
print(summary_stats)
