
################### Muskox analysis 2022 ######################################
################ Read in the data and load the packages #######################
rm(list=ls()) # clear the environment
setwd("C:/Users/elean/OneDrive - University of Calgary/Kutz Postdoc/Muskox body condition/Data analysis") # set the working directory

# load libraries
library(dplyr) # for data manipulation
library(devtools) # for model support
library(plspm)
library(ggplot2)
library(ggthemes)
library(gridExtra)


########### Read and handle the data ##########################################
# read the data
data <- read.csv("Muskox_master.csv")


# investigate the structure of the data
head(data) # shows the top 6 rows
names(data) # show the names of the columns

str(data)
data$sex <- as.factor(data$sex) # make the column a factor ie. categorical
data$age <- as.numeric(data$age) # make the column numeric ie. continuous
data$fetsex <- as.factor(data$fetsex) # make the column a factor ie. categorical

data$TotalTel <- as.numeric(data$Male_tel + data$Female_tel)
data$TotalMarsh <- as.numeric(data$Male_marsh + data$Female_marsh)

#  Omit rows without abomasum counts
data <- data[!is.na(data$TotalTel),]
#  Omit rows without abomasum counts
data <- data[!is.na(data$llarvae),]

median(data$TotalMarsh)
median(data$TotalTel)

# Add parasite species richness
data <-  data %>% mutate(count=rowSums(.[c(22,24,26:28)]!=0))

data %>% 
  group_by(sex) %>% summarise(se = sd(KFI))
sd(data$KFI)



ggplot(data, aes(y= KFI, x = count))+
  geom_point()+geom_smooth(method="lm")

ggplot(data, aes(y= TotalTel, x = sex, group = preg))+
  geom_boxplot()

ggplot(data, aes(y= KFI, x = age))+
  geom_point()+geom_smooth(method="lm") +geom_smooth(colour = "red")


########### PLS-PM All individuals ############################################

head(data)
names(data)
data$age <- as.numeric(data$age) 


all <- data[,c(2,3,4,8,27,28,18,19,22,23,24,25,26,5,6)]
head(all)
all <- na.omit(all)

## Subset for males
males <- subset(all, sex == "m")

#Create the inner matrix. ie. the model paths 
# rows of the inner model matrix: these are the latent variables (LVs)
Demographics <- c(0, 0, 0)
Parasites <- c(1, 0, 0)
Body_condition <- c(1, 1, 0)

# create path matrix and add column names
muskox_path <- rbind(Demographics, Parasites, Body_condition) # put the rows together
colnames(muskox_path) <-  rownames(muskox_path) # name the columns

muskox_path # check what has been made - this shows the matrix structure
par(mfrow=c(1,1))
innerplot(muskox_path)# Create the outer model list 
## This defines the columns which will be used to inform each aspect of the inner model 
names(males)
muskox_blocks <- list(3, c(5:12), 4) # outer mode, manifest variables (MVs)

# set the muskox modes, to assume measurement of latent 
# variables is in reflective mode
muskox_modes <- c("A", "A", "A")

### Run the pls_pm 
muskox_pls <- plspm(males, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)
plot(muskox_pls,cex.txt = 2) # plot the latent variable relationships 
plot(muskox_pls, what = "loadings")

## remove vars
muskox_blocks <- list(3, c(5,6,11,12), 4) # outer mode, manifest variables (MVs)

### Run the pls_pm 
muskox_pls <- plspm(males, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)
plot(muskox_pls,cex.txt = 2)

## Subset for females
females <- subset(all, sex == "f")


muskox_blocks <- list(3, c(5:12), 4) # outer mode, manifest variables (MVs)

# set the muskox modes, to assume measurement of latent 
# variables is in reflective mode
muskox_modes <- c("A", "A", "A")

### Run the pls_pm 
muskox_pls <- plspm(females, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)


#### remove vars
muskox_blocks <- list(3, c(5,6,11,12), 4) # outer mode, manifest variables (MVs)

### Run the pls_pm 
muskox_pls <- plspm(females, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)
plot(muskox_pls,cex.txt = 2)


##################### GROUP COMAPRISON #########################################
# rows of the inner model matrix: these are the latent variables (LVs)
Demographics <- c(0, 0, 0)
Parasites <- c(1, 0, 0)
Body_condition <- c(1, 1, 0)

# create path matrix and add column names
muskox_modes <- c("A", "A", "A")
muskox_path <- rbind(Demographics, Parasites, Body_condition) # put the rows togethermuskox_blocks <- list(3, c(5:12), 4) 
muskox_blocks <- list(3, c(5:12), 4) # outer mode, manifest variables (MVs)
muskox_pls <- plspm(all, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)

muskox_blocks <- list(3, c(5,6,11,12), 4) # outer mode, manifest variables (MVs)
muskox_pls <- plspm(all, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
gpa_boot = plspm.groups(muskox_pls, all$sex)
summary(muskox_pls)
muskox_pls$unidim
gpa_boot

summary(gpa_boot)
gpa_boot$global
gpa_boot$group1
gpa_boot$group2
plot(gpa_boot)



########### Females individuals; effect of reproduction #######################
females <- subset(females, age >= 3)
females <- females %>% 
  mutate(lactation = ifelse(lact == "y", 1,0),
         pregnant = ifelse(preg == "y", 1,0))
head(females)
str(females)

ggplot(females, aes(y=pregnant, x=TotalTel))+
  geom_point()+geom_smooth(method="glm", color="green", se=FALSE, 
                           method.args = list(family=binomial))

#Create the inner matrix. ie. the model paths 
# rows of the inner model matrix: these are the latent variables (LVs)
Demographics <- c(0, 0, 0, 0, 0)
Lactation <- c(1, 0, 0, 0, 0)
Parasites <- c(1, 1, 0, 0, 0)
Body_condition <- c(0, 1, 1, 0, 0)
Reproduction <- c(1, 0, 1, 1, 0)

muskox_path <- rbind(Demographics, Lactation, Parasites, Body_condition, Reproduction) # put the rows together
colnames(muskox_path) <-  rownames(muskox_path) # name the columns

muskox_path # check what has been made - this shows the matrix structure
par(mfrow=c(1,1))
innerplot(muskox_path)# Create the outer model list 

muskox_modes <- c("A", "A", "A", "A", "A")

names(females)
muskox_blocks <- list(3, 16, c(5:12), 4, 17) # outer mode, manifest variables (MVs)

muskox_pls <- plspm(females, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)
plot(muskox_pls,cex.txt = 2) # plot the latent variable relationships 
plot(muskox_pls, what = "loadings")

## remove vars
muskox_blocks <- list(3, 16, c(5,6,11,12), 4, 17) # outer mode, manifest variables (MVs)

muskox_pls <- plspm(females, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls)

muskox_pls$unidim
Paths = muskox_pls$path_coefs
arrow_lwd = 20 * round(abs(Paths), 2)
plot(muskox_pls,cex.txt = 2, arr.lwd = arrow_lwd) 

############ Fetus development ################################################
all <- data[,c(2,3,4,8,27,28,18,19,22,22,23,24,25,26,5,6,9,10,11,12,29)]
females <- subset(all, sex == "f")
females <- subset(females, age >= 3)
head(females)
pregs <- subset(females, preg == "y")
head(pregs)
pregs <- na.omit(pregs)

#
Age <- c(0, 0, 0, 0)
Parasites <- c(1, 0, 0, 0)
Body_condition <- c(0, 1, 0, 0)
Fetus_devel <- c(1, 1, 1, 0)

muskox_path <- rbind(Age, Parasites, Body_condition, Fetus_devel) # put the rows together
colnames(muskox_path) <-  rownames(muskox_path) # name the columns

muskox_path # check what has been made - this shows the matrix structure
par(mfrow=c(1,1))
innerplot(muskox_path)# Create the outer model list 

muskox_modes <- c("A", "A", "A", "A")

### remove vars
names(pregs)
str(pregs)
muskox_blocks <- list(3, c(5:6), 4, c(20,21)) # outer mode, manifest variables (MVs)

muskox_pls3 <- plspm(pregs, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls3)
plot(muskox_pls3,cex.txt = 2)
muskox_pls3$unidim

##### Fetus Sex
Age <- c(0, 0, 0, 0)
Parasites <- c(1, 0, 0, 0)
Body_condition <- c(0, 1, 0, 0)
Fetus_devel <- c(1, 1, 1, 0)
names(pregs)

muskox_modes <- c("A", "A", "A","A")
muskox_path <- rbind(Age, Parasites, Body_condition, Fetus_devel) # put the rows together
colnames(muskox_path) <-  rownames(muskox_path) # name the columns

muskox_path # check what has been made - this shows the matrix structure
par(mfrow=c(1,1))
innerplot(muskox_path)# Create the outer model list 

muskox_blocks <- list(3, c(5:6), 4, 17) # outer mode, manifest variables (MVs)

muskox_pls3 <- plspm(pregs, muskox_path, muskox_blocks, modes = muskox_modes, boot.val = TRUE, br = 200)
summary(muskox_pls3)

