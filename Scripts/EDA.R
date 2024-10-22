### Competing risks data generated from UNOS; there are 2 events (delta = 1, 2) and censoring (delta = 0)
### Cleaned/Prepared by Edem Defor/Dipankar Bandyopadhyay (Virginia Commonwealth Univ.) 
### Date: 09/27/2023


rm(list=ls())

install.packages("librarian")

librarian::shelf(cmprsk, dplyr, survminer, purrr, tidyr, devtools) 

 
# loading Libraries

library(cmprsk)     # fits F&G (Fine&Gray approach), but not suitable for handling large dataset

# This loads the version 1.1.2 of fastcmprsk package, which runs very fast! A really cool recent contribution 
devtools::install_github("erickawaguchi/fastcmprsk")


packageVersion("fastcmprsk") ## Need to be version 1.1.2
library(dplyr)
library(survival)
library(survminer)
library(purrr)
library(tidyr)


# setting working directory; use your own directory!

#setwd("I:/My Drive/OneDrive/Work/VCU_PhD/EdemDefor/EdemWork/EricaMoodie/")
setwd("C:/Users/misha/Documents/Thesis/kidney_data")

# read the dataset directly from the folder
df01 <- read.csv("cmprsk_2001.csv")
dim(df01)


str(df01)



######## Data cleaning ####

# identifying columns with more than 2 but less than 10 categories

      cols1 = mapply(table,df01) # table all the data
     
      index = NULL  # set to assign columns are >2 but <10 categories
      for(i in 1:ncol(df01)){
        if(length(cols1[[i]])>2 & length(cols1[[i]])<10){
          index = append(index, i)
        }
      }
      # convert columns with >2 but <10 columns to factor
      df2 =  df01
      df2[,index] <- lapply(df2[,index], factor)
      
      str(df2[,index]) # check if those columns are needed to be factors
      
      df2 = df2 %>%  mutate_at(vars(matches("delta")),as.character) # reconvert delta to numeric
      
      df2 = df2 %>% mutate_at(vars(matches("delta")),as.integer)  # reconvert delta to numeric
      
      str(df2) # check the data structure
      
# missing values      
df.complete <- df2[complete.cases(df2[, c("timeto","delta")]),]


## Please note that henceforth, we will be working with the "df.complete" dataset

## Now, we will generate proportion of outcomes, censoring and competing event
## The code below will give you the # censored (delta = 0), # (cause = 1), and # (Cause = 2)


df.complete %>%
  group_by(delta) %>%
  summarise(n=n())%>%
  mutate(freq = round((n / sum(n)),4))



## The code block below is to handle en bloc transplants (dual kidney); please ignore! 
df.enbloc = df.complete[df.complete$tx_enbloc==1,]
df.enbloc %>%
  group_by(delta) %>%
  summarise(n=n())%>%
  mutate(freq = round((n / sum(n)),4))
  
  
## Now, I am trying to compute the 2x2 table wrt, Donor HCV (DonorHCV) and Recipient HCV (RecHCV) status! HCV is Hep C virus
##Clearly, this should represent a very asymmetric table (no wonder), given that a lot of transplants 
##are expected to be from the HCV -ve donors to the HCV -ve recipients. In recent times, HCV+ donors are giving kidneys to normal/HCV- 
## recipients. This direction of research is continuing, but I do not find any nice stochastic formulation to this, if any,
## from a DTR, or optimization standpoint. This can likely be the starting point of our work. Furthermore, how to address the 
## clustering issue, i.e., subjects within centers? 


table(df.complete$DonorHCV, df.complete$RecHCV)



## Read some docs below to understand the HCV +/- donor story: 
## https://www.kidneyfund.org/all-about-kidneys/other-kidney-problems/hepatitis-c/hepatitis-c-and-kidney-transplants#
## https://pubmed.ncbi.nlm.nih.gov/33259125/
## https://jamanetwork.com/journals/jama/fullarticle/2795744



##### Now, we like to work on some simple CSH regressions; here we are just using covariates we like, but one can do this after testing the model #####
##### Covariate dictionary is also attached ########
##### Also note, we may NOT want to do CSH, and using Cox. We want to use AFT, and use CIFs (the Fine&Gray approach) 

## Cause 1 - graft failure

fit1 <- coxph(Surv(timeto, ifelse(delta == 1, 1, 0)) ~ donorage+DONORAA+genderd_num+hgt_don+wgt_don+cardarrest+
                +DONORHTN+DONORDM+DONCIG+DonorHCV+DONOR_SER+immuno_group
                +frailty(ctrcode),data=df.complete) 

summary(fit1)











