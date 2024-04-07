## Libraries
library(rjags)
library(renv)

setwd("/Users/johnnymosby/repos/thesis_demography/bayesianAnalysis/script")

data = read.csv("microdata.csv") |> 
    subset(select = c("m_age_at_birth", "f_age_at_birth",
                      "m_annual_salary_tm2", "f_annual_salary_tm2",
                      "gr_m_family_status",
                      "edu4_comp")) |> 
    na.omit()



modelAsText = "
model {
  
}
"
