# brms formula
list_formulas <- paste(paste("cbind(", paste(colnames(Y[-1]), collapse = " , "), ") ~ ", sep = ""), 
                       paste(colnames(X$x[,1:nComp]), collapse = " + "),  
                       '+ (1 | p | gr(site_ID, "normal")')#, "+ (1 | b | species_ID)")
paste(paste("cbind(", paste(colnames(Y[2]), collapse = " , "), ") ~ ", sep = ""), 
      paste(colnames(X$x[,1:nComp]), collapse = " + "),  
      "+ (1 | p | gr(site_ID))")

list_formulas <-bf(paste(colnames(Y[2]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "gaussian") +
  bf(paste(colnames(Y[3]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[4]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[5]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "gaussian") +
  bf(paste(colnames(Y[6]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[7]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[8]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[9]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "gaussian") +
  bf(paste(colnames(Y[10]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[11]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[12]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "gaussian") +
  bf(paste(colnames(Y[13]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[14]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[15]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[16]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student") +
  bf(paste(colnames(Y[17]), " ~ ", paste(colnames(X$x[,1:nComp]), collapse = " + "),"+ (1 | p | gr(site_ID))"), family = "student")
  


list_formulas <- paste(paste("cbind(", paste(colnames(Y[c(5,16)]), collapse = " , "), ") ~ ", sep = ""), 
                       paste(colnames(X$x[,1:nComp]), collapse = " + "),  
                       "+ (1 | p | gr(site_ID))")#, "+ (1 | b | species_ID)")

str.md <- make_stancode(list_formulas, data = traindata,family = student()) #, autocor = cor_arr())

data(oldcol, package = "spdep")


