
# NeuroVault 2447
# Contrast: noun verb judgement > baseline
library(httr)
library(jsonlite)
library(RNifti)
library(ARIbrain)

collection_id <- 2447
out_dir <- paste0("neurovault_", collection_id)
dir.create(out_dir, showWarnings = FALSE)

api_url <- paste0("https://neurovault.org/api/collections/", collection_id, "/images/")
res <- GET(api_url)
data <- fromJSON(content(res, "text", encoding="UTF-8"))


sel <- data$results$contrast_definition == "noun/verb judgement > baseline"
sel[is.na(sel)] <- FALSE
image_urls <- data$results$file[sel]
image_ids  <- data$results$id[sel]

n <- length(image_urls)
sids <- c(paste0("0",1:9),10:n)

for (i in seq(n)) {
  fname <- file.path(out_dir, paste0("tmap_", sids[i], ".nii.gz"))
  download.file(image_urls[i], destfile=fname, mode="wb")
}



df <- data$results[sel, ]
names(df)                 
table(df$analysis_level)

# We do not have a MASK, so we create one that filters only voxels with at least one non-zero statistic

n <- 10
sids <- c(paste0("0",1:9),10:n)

any_nonzero  <- NULL
any_infinite <- NULL

for (i in seq(n)) {
  file_path <- file.path(out_dir, paste0("tmap_", sids[i], ".nii.gz"))
  nii <- readNifti(file_path)
  arr <- as.array(nii)
  
  nonzero  <- (arr != 0)
  infinite <- is.infinite(arr)
  
  if (is.null(any_nonzero)) {
    any_nonzero  <- nonzero
    any_infinite <- infinite
  } else {
    any_nonzero  <- any_nonzero  | nonzero
    any_infinite <- any_infinite | infinite
  }
}

#The resulting mask_sum is a comprehensive union of all non-zero voxels found across the entire collection of t-maps. 
#This is often used as a way to define a minimal analysis region that includes all relevant activity found in the dataset.

mask_sum <- any_nonzero & !any_infinite

mask_nii <- asNifti(mask_sum*1, reference=nii)
writeNifti(mask_nii, file.path(out_dir, "mask.nii.gz"))

mask <- RNifti::readNifti("~/neurovault_2447/mask.nii.gz")
mask_v <- as.vector(mask)


# -----------------------------------------------
# P-VALUES FOR VOXELS INSIDE MASK 
# -----------------------------------------------

n <- 10
m <- sum(mask)

# Create matrix of t statistics, with n rows = subjects and m columns = voxels

copes <- matrix(0, nrow=n, ncol=m)
sids <- c(paste0("0",1:9),10:n)

for (i in seq(n)){
  fname <- paste0("~/neurovault_2447/tmap_", sids[i], ".nii.gz")
  copes[i,] <- as.vector(RNifti::readNifti(fname))[mask_v==1]
}

# Transform it into vector and compute p-values
# Note: the matrix is converted column-wise, so that the first x[1:n] correspond to voxel 1,
# x[(1:n) + n] to voxel 2, and in general
# x[(1:n) + (j-1)n] to voxel j

x <- as.vector(copes)

# Compute p-values approximating the t distribution with a standard normal
p <- 2*pnorm(abs(x), lower.tail=FALSE)
save(p, mask, file="words_pvals.RData")

# ----------------pmap for each subject-----------------#

for (i in seq(n)) {
  
  tfile <- file.path(out_dir, paste0("tmap_", sids[i], ".nii.gz"))
  tmap  <- readNifti(tfile)
  tarr  <- as.array(tmap)
  
  # two-tailed p-values
  parr <- 2 * (pnorm(abs(tarr), lower.tail = F))
  
  # optional: apply your mask
  parr[!mask_sum] <- NA
  
  pnii <- asNifti(parr, reference = tmap)
  
  writeNifti(
    pnii,
    file.path(out_dir, paste0("pmap_", sids[i], ".nii.gz"))
  )
}

# ------------------------------------------
# VISUALIZE
# ------------------------------------------

## ------------------------------------------
## EXAMPLE OF SUBJECT j 
## ------------------------------------------



sum(mask)
mask <- RNifti::readNifti("~/neurovault_2447/mask.nii.gz")
mask_v <- as.vector(mask)

subj1 <- RNifti::readNifti("~/neurovault_2447/tmap_01.nii.gz")
subj_v <- as.vector(subj1)[mask_v==1]

max(subj_v)

summary(subj_v)

hist(subj_v)
hist(subj_v[abs(subj_v)<=5])

p_subj1 <- 2*pnorm(abs(subj_v), lower.tail=FALSE)
hist(p_subj1)

subj1_p <- readNifti("~/neurovault_2447/pmap_01.nii.gz")
subj1_p_vec<- as.vector(subj1_p)[mask_v==1]

#sbj1 t map
vis_brain(subj1, mask, x=31, y=6, z=20)

#sbj 1 p map
vis_brain(subj1_p, mask, x=31, y=6, z=20, type = "abs")


############### Standard ARI ##################



# To apply standard ARI, we first compute test statistics and p-values
# for each voxel inside the mask

tmp <- pARI:::oneSamplePar(t(copes), alternative = "two.sided")
ts <- tmp$Test
p <- tmp$pv


# Then we pass from a vector X (either t statistics or p-values) to a map

from_v_to_map <- function(X, pvalue=FALSE){
  out <- mask_v
  out[mask_v==1] <- X
  
  tau <- ifelse(pvalue,1,0)
  out[mask_v==0] <- tau
  
  out <- array(out, dim=dim(mask))
  return(out)
}


tmap <- from_v_to_map(ts)
pmap <- from_v_to_map(p)


# Then we get super-threshold clusters, for t statistics over a certain threshold
# For simplicity, we select only clusters with a minimum size

get_cl <- function(thr=3.2, min_size=20){
  
  out <- ARIbrain::cluster_threshold(abs(tmap) > thr)
  
  tb <- table(out)
  to_exclude <- as.numeric(names(tb[tb < min_size]))
  
  ix <- out %in% to_exclude
  ix <- array(ix, dim=dim(out))
  out[ix] <- 0
  
  return(out)
}



cl1 <- get_cl(3.2)
table(cl1)

res <- ARIbrain::ARI(pmap, clusters=cl1, mask=mask)
res

# Remove row for cluster zero (i.e., voxels outside cluster)
res <- res[-nrow(res),]
res



# For the plot, we just assign the value of the TDP (ActiveProp)
# to each cluster

resmap <- array(0, dim=dim(cl1))


x <- rownames(res)
x <- as.numeric(sub("^cl", "", x))

for(i in seq_along(x)){
  resmap[cl1==x[i]] <- res[i,4]
}


resmap[mask==0] <- NA


visualizeBrain(resmap, x=31, y=6, z=20)

maxTDP <- max(resmap, na.rm = TRUE)
vis_brain2(resmap, mask, x=31, y=6, z=20, type="tdp", limits=c(0, maxTDP))

