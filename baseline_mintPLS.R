#test on mintPLS
y_obs <- augmented_matrix[!(colnames(augmented_matrix) %in% 
                              c("flightpath", "band_site", "band_species", paste("band_", seq(1,369), sep="")))] %>%
  unique
train_data <- augmented_matrix
X <- hiper_features(train_data, normalization = "norm2") #%>% #, normalization = "no") %>%
  #prcomp 
Y <- inner_join(train_data["individualID"], y_obs) #train_data[!grepl("band", names(train_data))] %>%
traindata <- data.frame(Y, X,#$x[,1:nComp], 
                        train_data$band_species, train_data$band_site) %>%
  data.frame %>%
  group_by(individualID) %>%
  filter(row_number()==ceiling(n()/2))
n_features <- dim(traindata)[2]-1
colnames(traindata)[(n_features):(n_features+1)] <- c("species_ID", "site_ID")

#plot reflectances
plot_data <- traindata %>% 
  dplyr::select(-one_of(c(colnames(Y), "site_ID", "species_ID"))) 
plot_data <- plot_data[-1] %>%
  t %>%
  data.frame
colnames(plot_data) = unlist(traindata[1]) # the first row will be the header
plot_data <- data.frame(bnd = 1:dim(plot_data)[1], plot_data)
ggdat <- tidyr::gather(plot_data, treeID, Reflectance,-bnd)

ggplot(ggdat[1:3530,], aes(x = bnd, y = Reflectance)) + 
  geom_line(aes(color = factor(treeID), alpha= 1), size = 0.2) +
  theme(legend.position="none")
  

library(mice)
imp <- mice(traindata, m = 1, print = FALSE)
dat <- complete(imp, 3)
mintPLS <- mixOmics::mint.spls(dat[colnames(dat) %in% colnames(X)], 
                               dat[colnames(dat) %in% colnames(Y[-1])], ncomp = 40, mode = "regression", 
                               dat$site_ID, scale = TRUE, tol = 1e-06, max.iter = 100,
                               near.zero.var = TRUE, all.outputs = TRUE)


test_mat <- readr::read_csv('Model_builder/dat/Crown_outBag.csv') 
test_mat$band_species <- factor(test_mat$band_species, levels = levels(augmented_matrix$band_species))
test_mat$band_site <- factor(test_mat$band_site, levels = levels(augmented_matrix$band_site))
test_mat <- test_mat #%>%
#dplyr::select(-one_of(c("foliarIronConc", "foliarPotassiumConc")))

Y_tst <- test_mat[!(colnames(test_mat) %in% 
                      c("flightpath", "band_site", "band_species", paste("band_", seq(1,369), sep="")))] 
y_test <- Y_tst %>% unique
X_tst <- hiper_features(test_mat, normalization = "norm2") #, normalization = "no")
pred <- predict(mintPLS, newdata = X_tst, test_mat$band_site)
ncomp = 10
mypls <- sapply(1:ncol(Y_tst[-1]), function(i) get_mod_r2(pred =  pred$predict[,,ncomp][,i], obs = as.matrix(Y_tst[,i+1])))
mypls

myR2 <- list()
Y_selected <- y_test #%>%
for(i in 1:ncol(Y_tst[-1])){
  y_hat <-  pred$predict[,,ncomp][,i] #*sd_traits[i] + mean_traits[i]
  
  bag_median <- data.frame(Y_tst[1], y_hat) %>%
    group_by(individualID) %>%
    summarize_all(mean)
  bag_median <- data.frame(Y_selected[1], Y_selected[,i+1]) %>%#*sd_traits[i] + mean_traits[i]) %>%
    inner_join(bag_median, by = "individualID")
  bag_median = bag_median[complete.cases(bag_median), ]
  myR2[[i]] <- get_mod_r2(bag_median[,2], bag_median[,3])
}
colnames(y_test[-1])
unlist(myR2)
