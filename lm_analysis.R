# Cacao Landmark Analysis for Wolcott et al. 2023, APPS
# Last updated on Apr 15 by K Wolcott

# If on Mac, first install XQuartz 2.7.11
## https://www.xquartz.org/releases/XQuartz-2.7.11.html

# Open terminal, cd to lm directory (one folder per specimen with multiple lm files), and remove lms only used for measurements, not gmm
# LM names to remove: OvaryH, StigLen, StigToStam, StyleH, ovary_style_stam, ovary_perimeter
#find . -name '*Ovary*' #-delete
#find . -name '*ovary*' #-delete
#find . -name '*Stig*' #-delete
#find . -name '*Style*' #-delete

# Load packages using pacman
#install.packages("pacman")
pacman::p_load(devtools, geomorph, rgl, purrr, dplyr, glue, abind, stringr,
               SlicerMorphR, ggplot2, ggforce, concaveman, Rcpp) 

# Set working directory  (define your data_wd below)
data_wd = "/Users/katherinewolcott/Documents/r/cacao"
lms_wd = paste(data_wd, "/lms", sep="")
setwd(lms_wd)

# Build df of all landmark files for each specimen
## Get paths for all specimens (1 folder each containing several landmark files)
dirs = list.dirs(".", full.names=FALSE, recursive=FALSE)
## Populate df with all landmark files per specimen
all_lms = lapply(dirs, function(y) {
    glue("Fetching landmark files for specimen: {y}")  
    path = paste(y, '/lms/', sep='')
    files <- dir(path, pattern = "*.json")
    glue("Compiling landmark files: {files}") 
    lms <- files %>% map_df(~data.frame(read.markups.json(
      file.path(path, .))) %>% mutate_if(is.numeric, as.character))
    rownames(lms) = NULL
    lms["spec_id"] = y
    rbind(lms)
})
## Optional: Check df to make sure no specimens are missing any lms
#for(i in 1:length(all_lms)) {        
  #print(i)
  #print(nrow(all_lms[[i]][1]))             
#}
## Export to file
df = bind_rows(all_lms) # Normalize dataframe with landmarks as rows, coords, spec_id, and lm_id as cols
write.table(df, "all_lms_wout-names.tsv", sep = "\t", row.names = FALSE)

# Add specimen and landmark names to all_lms dataframe
## Make list of specimen names
fpath_idx = lapply(dirs, function(y) {
    path = paste(y, '/lms/', sep='')
    files <- dir(path, pattern = "*.json")
    fpath_idx = list(path)
    })
fpath_idx = unlist(fpath_idx)
write.table(fpath_idx, "specimen_idx.tsv", sep = "\t")
## Add specimen and landmark names to all_lms dataframe
master_lm_df = read.table("all_lms_wout-names.tsv", sep = "\t", header = TRUE)
## Make list of landmark names (lm curve name, number of lm points)
## hacky step, but avoids file name irregularities resulting from duplicate naming in 3D slicer
lm_names = c(rep("fil_bot", 15), rep("fil_top", 10), rep("ligule_length", 20),
             rep("ligule_width", 15), rep("petal_curveL", 20), rep("petal_curveR", 20),
             rep("sepal_length", 20), rep("sepal_width", 10), rep("staminode_height", 5))
master_lm_df["lm_name"] = rep(lm_names, length(unique(df$spec_id)))
## Add numeric labels for 135 landmarks for each specimen 
master_lm_df["lm_id"] = rep(c(1:135), length(all_lms))
master_lm_df = master_lm_df[, c(4, 5, 6, 1, 2, 3)] # re-order columns for easier reading
## Optional: Check if there are any NA's
#master_lm_df[!complete.cases(master_lm_df),]
## Write to file File S2 - all_lms
write.table(master_lm_df, "all_lms.tsv", sep = "\t", row.names = FALSE)

# Format landmark coordinates as an array for geomorph
df_lms = df[c('X', 'Y', 'Z')] # Take only coodinates from all_lms df
df_lms_numeric = mutate_all(df_lms, function(x) as.numeric(as.character(x))) # Force as numeric
m = data.matrix(df_lms_numeric) # Make into an array
a = arrayspecs(A=m, p=135, k=3, sep=NULL) # Define array shapes

# Inspect data
a[,,1] # Inspect one specimen
a[1,,] # Inspect one landmark for all specimens
a[,1,] # Inspect landmark coord 1 of 3 for all specimens

# Define sliding landmarks for curves of flower
setwd(data_wd)
sliders = read.table("lm_curves.tsv", header=TRUE)
sliders = sliders[4:6]

# Make specimen info df to sort lms by groups (location, month, self compat, etc)
## Load in specimen id mapping file with names, collection info
spec_ids = read.table("map_specimen_ids.tsv", header=T, sep="\t")
spec_ids$spec_id<-tolower(spec_ids$spec_id) 
## Load in specimen idx file to add specimen names to numbers for interpreting lm analysis
setwd(lms_wd)
spec_nums = read.table("specimen_idx.tsv", sep = "\t", row.names=1)
spec_nums["spec_num"] = row.names(spec_nums)
spec_nums["spec_id"] = spec_nums$x
spec_nums["spec_id"] = str_replace(spec_nums[,3], '/lms/', '')
spec_nums$spec_id<-tolower(spec_nums$spec_id) 
## Combine files
spec_info <- merge(spec_ids, spec_nums)
colnames(spec_info)[7] <- 'tree_id'
## Set variables of interest as factors
spec_info$col_month = as.factor(spec_info$col_month) # collection month
spec_info$selfCompatible = as.factor(spec_info$selfCompatible) # self compatible
spec_info$col_loc = as.factor(spec_info$col_loc) # collection location
spec_info$tree_id = as.factor(spec_info$tree_id) # tree id

# Run GPA to align and scale specimens use Procrustes Distance for sliding
Y.gpa <- gpagen(a, curves = sliders, ProcD = TRUE)
summary(Y.gpa)
plot(Y.gpa)

# Get mean shape coordinates (reference shape) after GPA alignment
ref = mshape(Y.gpa$coords, na.action = 3)
plot(ref)

# Compare a specimen to mean consensus shape
plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
plotRefToTarget(ref,Y.gpa$coords[,,1],method="TPS")

# Plot outliers
outliers = plotOutliers(Y.gpa$coords)

# Find the specimen with the closest values to the mean shape
meanspec = findMeanSpec(Y.gpa$coords)

# Compare an individual specimen to the mean shape using different plotting methods
target = Y.gpa$coords[,,1] # Normal specimen 
## Thin	plate	spline
plotRefToTarget(ref, target, method="TPS")
## Lollipops
plotRefToTarget(ref,	target, method="vector")
## Points
plotRefToTarget(ref,	target, method="points")

# PCA, all specimens
pca = gm.prcomp(Y.gpa$coords)
sum = summary(pca)
prop_var = sum$PC.summary[2,1:2]

# Plot PCA in ggplot
dtp = data.frame('tree_id'=spec_info$tree_id, # Take first 2 components of PCA
                 'PC1'=pca$x[,1], 'PC2'=pca$x[,2]) 
## Plot parameters
p = ggplot(data = dtp, aes(x = PC1, y = PC2, col = tree_id)) + 
          # Add polygons around each tree
          geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=tree_id))+
          geom_point() +  
          theme_minimal()
## Add labels
p + labs(x=glue("PC1: {round(prop_var[1], 4)*100}%"), 
         y=glue("PC2: {round(prop_var[2], 4)*100}%"))

# Make deformation meshes to manually add to PCA plots outside of R
## Compare min and max from PCA to mean consensus shape
ref<-mshape(Y.gpa$coords) # mean shape
## Component 1
plotRefToTarget(ref, pca$shapes$shapes.comp1$min) # min shape to mean
plotRefToTarget(ref, pca$shapes$shapes.comp1$max) # max shape to mean
## Component 2
plotRefToTarget(ref, pca$shapes$shapes.comp2$min) # min shape to mean
plotRefToTarget(ref, pca$shapes$shapes.comp2$max) # max shape to mean

# Procrustes Anova
## Make a geomorph dataframe
gdf = geomorph.data.frame(Y.gpa, tree=spec_info$tree_id, 
                          col_loc=spec_info$col_loc, col_month=spec_info$col_month,
                          self_compat=spec_info$selfCompatible)
attributes(gdf)
## Run anova
### Null fit model (no group effects)
fit = procD.lm(coords ~ log(Csize), data = gdf)
aov = anova(fit)
summary(aov)

### Common fit model (group effects of tree_id, col_month)
fit.common = procD.lm(coords ~ log(Csize) + tree + col_month, data = gdf, SS.type = "III")
aov.common = anova(fit.common)
summary(aov.common)

### Common fit model sc (group effects of tree_id, self_compat, col_month)
fit.common.sc = procD.lm(coords ~ log(Csize) + self_compat + tree + col_month, data = gdf, SS.type = "III")
aov.common.sc = anova(fit.common.sc)
summary(aov.common.sc)

### Common fit model sc (group effects of tree_id, self_compat, col_month)
fit.common.sc = procD.lm(coords ~ log(Csize) + self_compat + tree + col_month, data = gdf, SS.type = "III")
aov.common.sc = anova(fit.common.sc)
summary(aov.common.sc)

# Compare fit models
anova(fit, fit.common, fit.common.sc, SS.type = "III")

# Pairwise test of unique vs unique sc
PW = pairwise(fit = fit.common.sc, fit.null = fit, groups = gdf$self_compat, covariate=gdf$Csize)
summary(PW, confidence = 0.95)

PW = pairwise(fit = fit.common.sc, fit.null = fit, groups = interaction(gdf$self_compat, gdf$tree))
summary(PW, confidence = 0.95)

## Export results
capture.output(summary(aov),file="anova.txt")

# Generate descriptive statistics for Table 1 - flower measurements
setwd(data_wd)
fpath = "measurements_raw.tsv"
meas = read.table(fpath, sep="\t", header=TRUE)
## Inspect raw data
head(meas)
str(meas)
## Drop columns with qualifying variables, keep only measurements
meas_cols = na.omit(meas[,4:18])
## Calculate min, mean, med, max, and sd
meas_stats = data.frame(lapply(meas_cols, min), row.names = "min")
meas_stats = rbind(meas_stats, data.frame(lapply(meas_cols, mean), row.names = "mean"))
meas_stats = rbind(meas_stats, data.frame(lapply(meas_cols, median), row.names = "median"))
meas_stats = rbind(meas_stats, data.frame(lapply(meas_cols, max), row.names = "max"))
meas_stats = rbind(meas_stats, data.frame(lapply(meas_cols, sd), row.names = "sd"))
## Calculate rsd (mean/sd)
t = data.frame(t(meas_stats))
t["rsd"] = t[,5]/t[,2]
## Export to file - Table 1
write.table(t, "Table_1-measurements.tsv", sep="\t")
