# libraries
library(VIM) # handling missing values
library(Amelia) # plot missing values
library(tidyverse) 
library(ggmap) # visualization of virus
library(corrplot) # plot of correlation coefficients
library(MASS)
library(broom) # model diagnostics of logistic regression 
library(effects) # plot of cordinal logistic regression
library(RColorBrewer)
# load data set
data <- read.csv('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 1 - Enteric Viruses\\Data\\Data for enteric virus MSc.csv')

# EDA
## summarize data
data.temp <- read.csv('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 1 - Enteric Viruses\\Data\\Data for enteric virus MSc.csv', stringsAsFactors = FALSE)
### age
p <- ggplot(data.temp, aes(x=CentrallyCity, y=Age, fill=CentrallyCity)) + geom_boxplot()

### Year of enrolment
sum <- summary(data)
table(data$Year.of.enrollment)
### the number of patients in each CentrallyCity
data.temp <- as.tibble(data.temp) 
num.of.patients <- data.temp %>%
  select(CentrallyCity, Age) %>%
  group_by(CentrallyCity) %>%
  summarise(number = n())

num.of.patients %>%
  mutate(CentrallyCity = fct_reorder(CentrallyCity, number)) %>%
  ggplot(aes(x=CentrallyCity, y=number)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_text(aes(label = number),
            hjust = -0.1,
            color = "#f68060",
            size = 3) +
  coord_flip() +
  xlab("") +
  theme_bw()

### the visualization of the distribution of virus
# before using ggmap, we need to upload the API key of Google map: register_google()
centre <- c(mean(range(data.temp3$LONGITUDE)), mean(range(data$LATITUDE)))
viet <- get_map(centre, zoom = 6, maptype = "hybrid")
ggmap(viet)

ggmap(viet) +
  geom_point(data = data.temp3, aes(x = LONGITUDE, y = LATITUDE),colour = "red", size = data.temp3$is.common)+
  geom_point(data = data.temp3, aes(x = LONGITUDE, y = LATITUDE),colour = "blue", size = data.temp3$is.uncommon,alpha = .5)
qmplot(LONGITUDE, LATITUDE, data = data.temp3, maptype = "toner-lite", color = CentrallyCity)

common_virus <- data.temp3[,c(73,75,76,78,83,86)]
uncommon_virus <- data.temp3[,c(56:72,74,77,79:82,84:85,87)]
data.temp3$is.common <- 12 - rowSums(common_virus)
data.temp3$is.uncommon <- 52 - rowSums(uncommon_virus)

table.virus <- data.temp3 %>%
  group_by(CentrallyCity) %>%
  summarise(common.virus = sum(is.common),uncommon.virus = sum(is.uncommon))
table.virus <- as.matrix(table.virus[,-1])
rownames(table.virus) <- c("Dak Lak","Dak Nong","Dong Thap","Khanh Hoa","Quang Binh","Quang Tri","Thua Thien - Hue")
table.virus <- t(table.virus)

coul <- brewer.pal(3, "Pastel2")
barplot(table.virus, 
        col=coul , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        legend=rownames(table.virus),
        xlab="Region",
        main = 'The distribution of viruses in different regions')





### correlation analysis

#### coinfections: age
cor(data$Age, data$is_coinf) # -0.3039761

#### coinfections: Gender
data.temp2 <- as.tibble(data)
data.temp2 %>%
  select(Gender,is_coinf) %>%
  group_by(Gender) %>%
  summarise(mean = mean(is_coinf))

### Finding Redundant Attributes
#### removing all variables with the same value:
data['Betapapillomavirus'] <- NULL
data['Circovirus'] <- NULL
data['Lymphocryptovirus'] <- NULL
data['Alphapolyomavirus'] <- NULL
#### correlation coefficient:
numeric_cor <- data[c('Age','DOB_Year','Length.of.stay','ContactDiar','NumberDiarEpi',
                      'HaemoglobinResult','WhiteCellsResult','NeutrophilsResult','LymphocytesResult',
                      'EosinophilsResult','PlateletsResult','Temp')]
M <- cor(numeric_cor)
corrplot.mixed(M, tl.cex = 0.5)
corrplot(M, method = "number") # Display the correlation coefficient

# check the quality of data set
for(i in 56:87){
  data[,i][data[,i]==2]=0
}
rowSums(data[,56:87]) == data$is_coinf
## whether AdminDate is later than Date.of.hospital.entry
sum(as.Date(data$Date.of.hospital.entry, format='%d/%m/%Y') > 
  as.Date(data$AdminDate, format='%d/%m/%Y')) == 0
## whether Length.of.stay is equal to DateDischOrDeath minus Date of hospital entry
stay <- as.Date(data$DateDischOrDeath, format='%d/%m/%Y')-as.Date(data$Date.of.hospital.entry, format='%d/%m/%Y')
sum(stay != data$Length.of.stay) == 0
## whether KnownTemp is consistent to Temp
(is.na(data$KnownTemp) | data$KnownTemp==2) == is.na(data$Temp)


# Missing Values
# missmap: all variables
missmap(data, main = 'missmap of all variables')
# missmap: missed variables
miss.data <- data[c('is_coinf','NumberDiarEpi','HaemoglobinResult','WhiteCellsResult',
                  'NeutrophilsResult','LymphocytesResult','EosinophilsResult','PlateletsResult',
                  'KnownTemp','Temp','Ct.value', 'Ct.value.1', 'Ct.value.2', 'Ct.value.3', 
                  'Ct.value.4', 'Ct.value.5', 'Ct.value.6')]

missmap(miss.data, main = 'missmap of variables with NA')

# is_coinf
## missing data/how many varibales have missing values
missing1 <- colnames(data)[colSums(is.na(data)) > 0]
length(missing1)
## how many cases have missing values
missing2 <- rownames(data)[rowSums(is.na(data)) > 0]
length(missing2)
## impute the missing values
data$is_coinf[is.na(data$is_coinf)] = 0

# NumberDiarEpi
length(data$NumberDiarEpi[is.na(data$NumberDiarEpi)])
barplot(data$NumberDiarEpi, col = "red",border="red", main = 'NumberDiarEpi', xlab = 'Patients', ylab = 'Number of diarrhoel episodes')
## impute the missing values
data$NumberDiarEpi[is.na(data$NumberDiarEpi)]=mean(data$NumberDiarEpi,na.rm = T)
barplot(data$NumberDiarEpi, col = "blue",border="blue", main = 'NumberDiarEpi after imputation', xlab = 'Patients', ylab = 'Number of diarrhoel episodes')

# Blood test results
length(data$HaemoglobinResult[is.na(data$HaemoglobinResult)])
length(data$WhiteCellsResult[is.na(data$WhiteCellsResult)])
length(data$NeutrophilsResult[is.na(data$NeutrophilsResult)])
length(data$LymphocytesResult[is.na(data$LymphocytesResult)])
length(data$EosinophilsResult[is.na(data$EosinophilsResult)])
length(data$PlateletsResult[is.na(data$PlateletsResult)])

which(is.na(data$HaemoglobinResult))
which(is.na(data$WhiteCellsResult))
which(is.na(data$NeutrophilsResult))
which(is.na(data$LymphocytesResult))
which(is.na(data$EosinophilsResult))
which(is.na(data$PlateletsResult))

miss.data2 <- data[c('HaemoglobinResult','WhiteCellsResult','NeutrophilsResult',
                    'LymphocytesResult','EosinophilsResult','PlateletsResult')]
missmap(miss.data2, main = 'missmap of blood test results')
## remove missed case
data <- data[-704,]
## impute missing values
data$EosinophilsResult[is.na(data$EosinophilsResult)]=mean(data$EosinophilsResult,na.rm = T)
data$PlateletsResult[is.na(data$PlateletsResult)]=mean(data$PlateletsResult,na.rm = T)

# Temperature
length(data$Temp[is.na(data$Temp)])
length(data$KnownTemp[is.na(data$KnownTemp)])
## impute the missing values
data$KnownTemp[is.na(data$Temp)] = 2
data$Temp[is.na(data$Temp)]=mean(data$Temp,na.rm = T)

# Real-time PCR results
PCR <- c('Ct.value', 'Ct.value.1', 'Ct.value.2', 'Ct.value.3', 'Ct.value.4', 'Ct.value.5', 'Ct.value.6')
for(i in PCR){
  num <-  length(data[i][is.na(data[i])])
  cat(num, num/nrow(data),'\n')
}
## replace missing values with indicators
for(j in PCR){
  key1 <- is.na(data[j])
  key2 <- !is.na(data[j])
  data[j][key1] = 0
  data[j][key2] = 1
}

# ContactDiar
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
data$ContactDiar[data$ContactDiar == 9] <- getmode(data$ContactDiar)

# BloodStool
data$BloodStool[data$BloodStool == 9] <- getmode(data$BloodStool)

# MucoidStool
data$MucoidStool[data$MucoidStool == 9] <- getmode(data$MucoidStool)

# AbdominalPain
data$AbdominalPain[data$AbdominalPain == 9] <- getmode(data$AbdominalPain)

# ThreeDaysFever
data$ThreeDaysFever[data$ThreeDaysFever == 9] <- getmode(data$ThreeDaysFever)


# outliers
## Discover outliers with visualization tools
healthy <- data.temp3[is.na(data.temp3$is_coinf),]
unhealthy <- data.temp3[!is.na(data.temp3$is_coinf),]
OutVals1 <- boxplot(healthy$NumberDiarEpi)$out
which(healthy$NumberDiarEpi %in% OutVals1)
boxplot(healthy$Temp)$out

OutVals2 <- boxplot(unhealthy$NumberDiarEpi)$out
which(unhealthy$NumberDiarEpi %in% OutVals2)
boxplot(unhealthy$Temp)$out

# Creating New Variables in an existing dataset
## is.same.day
data$is.same.day <- data.temp$Date.of.hospital.entry == data.temp$AdminDate 
table(data$is.same.day)
## is.common.virus
virus.name <- names(data)[56:87]
name <- paste("is.common.", virus.name, sep="")
common.virus <- c('Rotavirus','Norovirus','Sapovirus','Kobuvirus',
                  'Mastadenovirus','Mamastrovirus')
common.name <- paste("is.common.",common.virus,sep = '')
for(i in name){
  if(i %in% common.name){
    data[i] <- 1
  }else{
    data[i] <- 0
  }
}

## Prevalence of viruses
barplot(colSums(data[,56:87])/nrow(data),main = 'Prevalence of viruses', xlab = 'Viruses', ylab = 'Rate')
colSums(data[,56:87])[which.max(colSums(data[,56:87]))]
colSums(data[,56:87])[which.min(colSums(data[,56:87]))]
## the frequency of different coinfections
data.temp2 <- read.csv('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 1 - Enteric Viruses\\Data\\Data for enteric virus MSc.csv')
data.temp2 <- data.temp2[,c(56:87)]
test <- data.temp2 %>%
  group_by_all() %>% summarise(COUNT = n())
test$row_sum = rowSums(test[,c(-33)])
## the prevalence of different common enteric viruses
data.temp3 <- read.csv('C:\\Users\\1ia0\\OneDrive - University of Edinburgh\\Dissertation (SDS)\\Session 1 - Enteric Viruses\\Data\\Data for enteric virus MSc.csv')

round(sum(data.temp3$Rotavirus == 1)/nrow(data.temp3),3)
round(sum(data.temp3$Norovirus == 1)/nrow(data.temp3),3)
round(sum(data.temp3$Sapovirus == 1)/nrow(data.temp3),3)
round(sum(data.temp3$Kobuvirus == 1)/nrow(data.temp3),3)
round(sum(data.temp3$Mastadenovirus == 1)/nrow(data.temp3),3)
round(sum(data.temp3$Mamastrovirus == 1)/nrow(data.temp3),3)



# modelling

## Poisson regression
plot(density(data$is_coinf))
hist(data$is_coinf)

### model fitting
data.ordinal <- data[,c(31:36,38:40,84)]

data.ordinal <- within(data.ordinal, {
  Tap <- factor(Tap)
  Well <- factor(Well)
  Rain <- factor(Rain)
  River <- factor(River)
  Pond <- factor(Pond)
  Bottled <- factor(Bottled)
  KeepAnimal <- factor(KeepAnimal, levels = c(2,1))
  KillingAnimal <- factor(KillingAnimal, levels = c(2,1))
  EatCookRawMeat <- factor(EatCookRawMeat, levels = c(2,1))
  is_coinf <- as.factor(is_coinf)
})

model.ordinal <- polr(formula = is_coinf ~ ., data = data.ordinal, Hess = TRUE)
sel.back.ordinal <- stepAIC(model.ordinal, direction = 'backward')

(ctable <- coef(summary(sel.back.ordinal)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))

### plots
plot(Effect(focal.predictors = "Well",sel.back.ordinal))
plot(Effect(focal.predictors = c("Well", "EatCookRawMeat"),sel.back.ordinal))

## whether/how coinfection impacts the disease severity
severity <- c('BloodStool', 'MucoidStool', 'NumberDiarEpi', 'AbdominalPain', 'ThreeDaysFever','Temp')
severity.df <- data[severity]
coinfection.df <- data[,56:83]
### logistic regression for discrete variables
severity.df$BloodStool[severity.df$BloodStool == 2] <- 0
BloodStool.df <- data.frame(y = as.factor(severity.df$BloodStool), x = data$is_coinf)
regr.BloodStool <- glm(y~., family = 'binomial', data = BloodStool.df)
summary(regr.BloodStool)

severity.df$MucoidStool[severity.df$MucoidStool == 2] <- 0
MucoidStool.df <- data.frame(y = as.factor(severity.df$MucoidStool), x = data$is_coinf)
regr.MucoidStool <- glm(y~., family = 'binomial', data = MucoidStool.df)
summary(regr.MucoidStool)

severity.df$AbdominalPain[severity.df$AbdominalPain == 2] <- 0
AbdominalPain.df <- data.frame(y = as.factor(severity.df$AbdominalPain), x = data$is_coinf)
regr.AbdominalPain <- glm(y~., family = 'binomial', data = AbdominalPain.df)
summary(regr.AbdominalPain)

severity.df$ThreeDaysFever[severity.df$ThreeDaysFever == 2] <- 0
ThreeDaysFever.df <- data.frame(y = as.factor(severity.df$ThreeDaysFever), x = data$is_coinf)
regr.ThreeDaysFever <- glm(y~.,family = 'binomial', data = ThreeDaysFever.df)
summary(regr.ThreeDaysFever)

### diagnostics of logistic regression model 

#### Influential values
Influential.check <- function(model){
  plot(model, which = 4, id.n = 3)
  temp.data <- augment(model) %>% 
    mutate(index = 1:n()) 
  p <- ggplot(temp.data, aes(index, .std.resid)) + 
    geom_point(aes(color = y), alpha = .5) +
    theme_bw()
  print(p)
  temp.data %>% 
    filter(abs(.std.resid) > 3)
}

Influential.check(regr.BloodStool)
Influential.check(regr.MucoidStool)
Influential.check(regr.AbdominalPain)
Influential.check(regr.ThreeDaysFever)

### linear regression model for continuous variables
NumberDiarEpi <- severity.df$NumberDiarEpi
NumberDiarEpi.df <- data.frame(y = NumberDiarEpi, x = data$is_coinf)
lm.NumberDiarEpi <- lm(y~.,data = NumberDiarEpi.df)
summary(lm.NumberDiarEpi)

Temp <- severity.df$Temp
Temp.df <- data.frame(y = Temp, x = data$is_coinf)
lm.Temp <- lm(y~.,data = Temp.df)
summary(lm.Temp)

## The analysis of risk factors associated with the common viruses
risk.df <- data[c(31:36,38:40)]

risk.common.virus <- function(x){
  virus <- data[x]
  virus.df <- cbind(virus,risk.df)
  formu <- paste(x,'~.',sep = '')
  regr <- glm(formu, family = 'binomial', data = virus.df)
  sel.back.virus <- stepAIC(regr, direction = 'backward',trace = 0)
  summary(sel.back.virus)
}

common.virus <- c('Rotavirus','Norovirus','Sapovirus','Kobuvirus','Mastadenovirus',
                  'Mamastrovirus')


risk.common.virus(common.virus[1])
risk.common.virus(common.virus[2])
risk.common.virus(common.virus[3])
risk.common.virus(common.virus[4])
risk.common.virus(common.virus[5])
risk.common.virus(common.virus[6])