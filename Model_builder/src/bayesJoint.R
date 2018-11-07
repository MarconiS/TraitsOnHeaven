library(rstan)
library(brms)
library(tidyverse)
coreNum <- 2

nComp = 20
augmented_matrix <- readr::read_csv('Model_builder/dat/Crown_norm.csv') 
set.seed(1)
augmented_matrix$band_species <- factor(augmented_matrix$band_species)
augmented_matrix$band_site <- factor(augmented_matrix$band_site)
augmented_matrix <- augmented_matrix #%>%
  #dplyr::select(-one_of(c("foliarIronConc", "foliarPotassiumConc")))
  #dplyr::select(c("foliarIronConc", "foliarPotassiumConc"))


source("Model_builder/src/utility.R")
y_obs <- augmented_matrix[!(colnames(augmented_matrix) %in% 
                              c("flightpath", "band_site", "band_species", paste("band_", seq(1,369), sep="")))] %>%
  unique
train_data <- augmented_matrix[complete.cases(augmented_matrix[colnames(augmented_matrix) %in% 
                                          paste("band_", seq(1,369), sep="")]), ]
X <- hiper_features(train_data, normalization = "norm2")  #%>% #, normalization = "no") %>%

#  prcomp 
Y <- inner_join(train_data["individualID"], y_obs) #train_data[!grepl("band", names(train_data))] %>%
#data.frame(Y[1], X) %>% dplyr::select(-one_of(colnames(Y[-1]))) %>% plot_spectra


traindata <- data.frame(Y, X, #$x[,1:nComp], 
                   train_data$band_species, train_data$band_site) %>%
  data.frame %>%
  group_by(individualID) %>%
  filter(row_number()==ceiling(n()/2))
n_features <- dim(traindata)[2]-1
colnames(traindata)[(n_features):(n_features+1)] <- c("species_ID", "site_ID")
traindata %>% dplyr::select(-one_of(colnames(y_obs[-1]))) %>% plot_spectra
pls_transform(traindata)
  
# i = 2
# ggplot(data = y_obs, aes_string(x = colnames(y_obs[i]))) + geom_histogram()
# i = i+1
# fill NAs
library(mice)
imp <- mice(traindata, m = 1, print = FALSE)
list_formulas <- paste(paste("cbind(", paste(colnames(Y[-1]), collapse = " , "), ") ~ ", sep = ""),
                       paste(colnames(X), collapse = " + "))#, "+ (1 | p | gr(site_ID))")

# list_formulas <- paste(paste("cbind(", paste(colnames(Y[-1]), collapse = " , "), ") ~ ", sep = ""), 
#                        #paste(colnames(X$x[,1:nComp]), collapse = " + "),  
#                        "(1 | p | gr(species_ID : site_ID))")#, "+ (1 | b | species_ID)")

fit <- brm_multiple(list_formulas, data = imp, cores =coreNum, chains = 2, iter = 200,  
                    set_prior("horseshoe()", class = "b"), 
                    #set_prior("lasso()", class = "cor"),
                    #set_rescor(TRUE), 
                    #autocor = cor_ar(~ 1 | site_ID),
                    #set_prior("lasso()", class = "ar"),
                    #set_prior("<prior>", class = "Intercept"),
                    #set_prior("<prior>", class = "sd", group = "<group>"),
                    family=brmsfamily("student"))

saveRDS(fit, "./sp_site_bjrm.rds")
sum(round(fit$rhats, 2) > 1.1)
#summary(fit)

#create test set
test_mat <- readr::read_csv('Model_builder/dat/Crown_outBag.csv') 
test_mat$band_species <- factor(test_mat$band_species, levels = levels(augmented_matrix$band_species))
test_mat$band_site <- factor(test_mat$band_site, levels = levels(augmented_matrix$band_site))
#test_mat <- test_mat #%>%
  #dplyr::select(-one_of(c("foliarIronConc", "foliarPotassiumConc")))

Y_tst <- test_mat[!(colnames(test_mat) %in% 
                              c("flightpath", "band_site", "band_species", paste("band_", seq(1,369), sep="")))] 
y_test <- Y_tst %>% unique
X_tst <- hiper_features(test_mat, normalization = "norm2") #, normalization = "no")
#X_tst <- predict(X, X_tst)
#Y_tst <- inner_join(test_mat["individualID"], y_test)  #train_data[!grepl("band", names(train_data))] %>%

test_data <- data.frame(Y_tst[1], X_tst, 
                   test_mat$band_species, test_mat$band_site) %>%
  data.frame 

n_features <- dim(test_data)[2]-1
colnames(test_data)[(n_features):(n_features+1)] <- c("species_ID", "site_ID")
posterior_dat <- predict(fit, newdata = test_data)

#Y_selected <- y_test #%>%
  #select(-one_of(c("foliarIronConc", "foliarPotassiumConc")))#[c(1,2,3,5,7, 8,9,15, 16, 18,20)]
#retrotransform
source("Model_builder/src/utility.R")
myR2 <- list()
for(i in 1:dim(posterior_dat)[3]){
  y_hat <-  (posterior_dat[,1,i]) #*sd_traits[i] + mean_traits[i]
  bag_median <- data.frame(Y_tst[1], Y_tst[(1+i)], y_hat) %>%
    group_by(individualID) %>%
    summarize_all(median)
  # bag_median <- data.frame(Y_selected[1], Y_selected[,i+1]) %>%#*sd_traits[i] + mean_traits[i]) %>%
  #   inner_join(bag_median, by = "individualID")
  bag_median = bag_median[complete.cases(bag_median), ]
  myR2[[i]] <- get_mod_r2(pred = unlist(bag_median[,2]), obs = unlist(bag_median[,3]))
}

myR2 <- data.frame(t(round(unlist(myR2), 2)))
colnames(myR2) <- colnames(y_test[-1])
myR2
# test_data <- data.frame(Y_tst, X_tst,#[,1:nComp], 
#                         test_mat$band_species, test_mat$band_site) %>%
#   data.frame 
# 
# n_features <- dim(test_data)[2]-1
# colnames(test_data)[(n_features):(n_features+1)] <- c("species_ID", "site_ID")
# R2 <- bayes_R2(current_bjrm_site, newdata = test_data[complete.cases(test_data), ])
# R2
readr::write_csv((myR2), "./myR2.csv")
#launch_shinystan(fit)

# plot(fit, pars = "nitrogenPercent")
# readr::write_csv(data.frame(R2), "./bayesR2.csv")

