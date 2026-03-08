
# ------------------------------------------
# COMPUTE TDP BOUNDS
# ------------------------------------------

library(RNifti)
require(hommel)

# I started with a matrix where n rows = subjects and m columns = voxels
# Then I trasformed the matrix into a vector

# For voxel i, the set jx of all subjects for that voxel is
# jx = (1:n) + (i-1)n


m <- sum(as.vector(mask)) # voxels inside the mask 
n <- length(p) / m # subjects

# Prepare the data
?hommel
hom <- hommel(p)

# For each voxel, track the number of active subjects (TD): "With 1-alpha confidence, there are at least $X$ subjects active in this voxel."
?hommel::discoveries
alpha <- 0.05
out <- rep(0,m)

for(j in 1:m){
  if(!(j %% 1000)){print(paste0(j,"/",m))}
  k <- (j-1)*n
  ix <- ((1:n)+k)
  out[j] <- hommel::discoveries(hom, ix=ix, alpha=alpha)
}

save(out, file="words_res.RData") #tdp's


# RECONSTRUCT THE NIFTI FILE

# Shiny App needs
consistency_img <- mask
consistency_img[] <- NA
# Fill the 1s in the mask with our 'out' values
consistency_img[mask == 1] <- out
writeNifti(consistency_img, "consistency_lower_bound_masked.nii")


# ------------------------------------------
# PLOT THE RESULTS
# ------------------------------------------

load("words_pvals.RData")
load("words_res.RData")

mask_v <- as.vector(mask)

# Transform results (out) into array, adding voxels outside mask as NA
res <- mask_v
res[mask_v > 0] <- out
res[mask_v == 0] <- NA
res <- array(res, dim=dim(mask))

dim(res)

vis_brain(res, mask, x=31, y=6, z=20, type = "abs")

table(as.vector(out))


# ------------------------------------------
# CLUSTERS OF VOXELS
# ------------------------------------------

# Get clusters of voxels that are active in at least gamma subjects (at least half of the sbj, at least 3,4 etc. try)

require(ARIbrain)
?ARIbrain::cluster_threshold

#gamma <- 5 #Only show me parts of the brain I am confident that at least 5 out of 10 subjects are active


make_cl <- function(gamma, min_size = 10){
  # 1. Create binary mask for voxels meeting the gamma criteria
  cond <- (res >= gamma)
  cond[is.na(cond)] <- FALSE
  
  if (sum(cond) < min_size) {
    return(res * 0) 
  }
  
  # 2. Identify clusters
  cl <- cluster_threshold(cond) 
  
  # 3. Filter clusters by size
  cl_sizes <- table(cl)
  # Keep only IDs where the count is >= min_size (excluding background 0)
  valid_clusters <- names(cl_sizes[cl_sizes >= min_size])
  valid_clusters <- valid_clusters[valid_clusters != "0"]
  
  # 4. Zero out small clusters
  cl_cleaned <- cl
  cl_cleaned[!cl %in% valid_clusters] <- 0
  
  # 5. Return the gamma value for significant clusters
  cl_gamma <- (cl_cleaned > 0) * gamma
  return(cl_gamma)
}



# Initialize with a size threshold of 10 voxels
cl_final <- make_cl(gamma = 1, min_size = 10)

n <- 10
for(gamma in 2:n){
  cl_current <- make_cl(gamma, min_size = 10)
  cl_final <- pmax(cl_final, cl_current, na.rm = TRUE)
}


visualizeBrain_01(cl_final, mask, x=31, y=6, z=20)


table(as.vector(cl_final))



# ------------------------------------------
# SUMMARY TABLE FOR HOMMEL RESULTS
# ------------------------------------------

# 'out' contains the TD bound for each voxel in the mask
# n is the number of subjects (usually 10)
n_subjects <- length(p) / m 

# 1. Count how many voxels have exactly s subjects active
# factor() ensures we see 0 to n, even if a count is zero
exact_counts <- table(factor(out, levels = 0:n_subjects))

# 2. Calculate "At Least" counts (Cumulative)
# We sum from the highest number of subjects down to zero
at_least_counts <- rev(cumsum(rev(as.vector(exact_counts))))

# 3. Calculate Percentages of the brain mask
percent_mask <- round((at_least_counts / length(out)) * 100, 2)

# Create the Data Frame
hommel_summary <- data.frame(
  s = 0:n_subjects,
  Voxels_Exact = as.vector(exact_counts),
  Voxels_At_Least = at_least_counts,
  Percent_of_Mask = percent_mask
)


print(hommel_summary)

