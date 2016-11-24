##first exploration meta-analysis#####

#clean workspace if needed###
rm(list = ls())

setwd("/Users/irishendriks/Documents/cursos/curso_UIB/MSc_Ecolog??a_marina/2016-17/meta-analysis Elena")

#load database
read.csv("Database.csv",header=TRUE,sep=",",dec=".")->data
head(data)
no.<-as.factor(data$no.)
N<-as.numeric(data$N)



#plot for visual inspection
plot(Sunscreen_ABB~Concentration.1, data=data)
Sunscreen_ABB

hist(data$Concentration.1)

# Simple plot to see if we have the same range of concentrations for all
boxplot(data$UNITS_ug_L~data$Sunscreen_ABB)
boxplot(data$UNITS_ug_L~data$group)
boxplot(data$UNITS_ug_L~data$Process)

####METAFOR####meta analysis
library(RCurl)
library(bitops)
library(metafor)
library(Formula)

data <- escalc(n1i = N, n2i = N, m1i = T_response, m2i = C_response, 
                  sd1i = T_var, sd2i = C_var, data = data, measure = "SMD", 
                  append = TRUE)
ma_model_1 <- rma(yi, vi, data = data)
summary(ma_model_1)
ma_model_2 <- rma(yi, vi, data = data, method="HE") #Hedges
summary(ma_model_2)

#forest plot
forest(ma_model_1, slab = paste(data$no., as.character(data$group), sep = ", "))

# Simple plot
boxplot(data$C_response~data$T_response, main="responses", xlab="Control response", ylab="Treatment response")


#logoddsratio
ES <- escalc(ai = T_response, bi = C_response, ci = T_var, di = C_var,
             data = data,
             measure = "OR")
cbind(ES$yi, ES$vi)

#random effect model, mean difference
result.md <- rma(m1 = T_response, m2 = C_response,
                 sd1 = sqrt(T_var), sd2 = sqrt(C_var),
                 n1 = N, n2 = N,
                 method = "REML", measure = "MD",
                 data = data)
summary(result.md)

#Get Between-Study Variance & Its Error
result.md$tau2
result.md$se.tau2
result.md$tau2 + 1.96 * c(-1, 1) * result.md$se.tau2 #95% CI


### meta-analysis of log ratio of means using a random-effects model

logratio_random <- rma(yi, vi, method="DL", data=data)
logratio_random
summary(logratio_random)
plot(logratio_random)

#In order to visualize the results you can create 
#a forest-plot using the forest() function.

forest(logratio_random, slab = paste(data$group), sep = ", ")

#### ahora para method "FE" FIXED EFFECTS 

rma_model_3 <- rma(yi, vi, method = "FE", data = data)
rma_model_3

###### Calcular el summary effect size por grupos . Method="DL" = DerSimonian-Laird estimator ####
### Bivalves
bivalves <- subset (data, group == "Bivalves")
bivalves_es <- escalc(n1i = N, n2i = N, m1i = T_response, m2i = C_response, 
               sd1i = T_var, sd2i = C_var, data = bivalves, measure = "SMD", 
               append = TRUE)
rma_model_4 <- rma(yi, vi, data = bivalves, method="DL")
summary(rma_model_4)
names(bivalves)
forest(rma_model_4, slab = paste(bivalves$group, sep = ", "))
title("Bivalves Random-effect model")

###### Microalgae, Fish and corals too few observations
algae <- subset (data, group == "Microalgae")
algae_es <- escalc(n1i = N, n2i = N, m1i = T_response, m2i = C_response, 
                      sd1i = T_var, sd2i = C_var, data = algae, measure = "SMD", 
                      append = TRUE)
rma_model_5 <- rma(yi, vi, data = algae, method="DL")
summary(rma_model_5)
forest(rma_model_5, slab = paste(algae$group, sep = ", "))
title("Algae Random-effect model")

# End ####
