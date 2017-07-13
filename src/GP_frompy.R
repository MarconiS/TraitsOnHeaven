library(rPython)

# send data frames to python
python.assign('test_spectra', test.data[,8:dim(test.data)[2]])
python.exec
python.assign('train_spectra', r_df1)
python.assign('test_traits', r_df1)
python.assign('train_traits', r_df1)
python.assign('specificTrait', r_df1)


# Load/run the main Python script
system(paste("python GetGP.py", args, sep = ""))

# Get the variable
new_subs_data <- python.get("new_subs")

# Load/run re-fecth script
python.load("RefreshNewSubs.py")

# Get the updated variable
new_subs_data <- python.get("new_subs")

head(new_subs_data)