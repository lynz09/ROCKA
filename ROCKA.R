# ROCKA
library(dplyr)
library(fpc)
library(mice)
library(ggplot2)
library(TTR)

setClass(Class = "Extract_baseline",
         representation(
           base_line = "data.frame",
           residual = "data.frame"
         )
)

setClass(Class = "ROCKA_cluster",
         representation(
           cluster_res = "numeric",
           centroid_dat = "data.frame"
         )
)

#########
# Preprocessing including missing values interpolation and standardization
# Baseline Extraction 
# including smoothing extreme values ,such as noises and anomalies
Smooth_anomaly <- function(dt){
  # Input:
  # dt: a dataframe with column id,days,y
  # method:a method to deal with missing data from R package 'mice'
  # remove the top 5% data which deviate the most from the mean value
  # and then use  interpolation methods to fill them
  # Output: a dataframe with column id,days,y
  if (anyNA(dt)==TRUE) {
    stop("Found NA in dataframe.")
  }
  if(!all(c('id','days','y') %in% names(dt))) {
    stop("Data have no columns 'id' or 'days' or 'y'. Please check your data.")
  }
  dt$days <- as.Date(dt$days)
  dt1 <-  dt %>% group_by(id) %>% 
    mutate(y_sd=scale(y),y_miu=abs(y-mean(y)))
  dt_n <- dt %>% group_by(id) %>% summarise(n=length(days))
  dt2 <- dt1 %>%  left_join(dt_n,by='id') %>% 
    arrange(id,desc(y_miu)) %>% 
    filter(row_number()<=round(n*0.05))
  dt3 <- dt1 %>% merge(dt2,by=c('id','days'),all.x = TRUE)
  dt3$y.x[!is.na(dt3$y.y)] <- NA # dataframe with removed anomalies
  wide <- dt3[,1:3] %>% setNames(nm=c('id','days','y')) %>% group_by(id) %>% 
    arrange(id,days) %>%
    dplyr::do({
      wide_ip <- mice::mice(data.frame(days=.$days,y=.$y),method = 'norm.predict')  # interpolation
      complete(wide_ip)
    })
  wide <- as.data.frame(wide)
  
 return(wide)
}
#########

#########
# Extract Baseline
Extract_Baseline <- function(dt,w=w){
  # separate curve into two parts:baseline and residuals
  # method:moving average with a sliding window  of length W
  # Input:dt ,a dataframe with column id,days,y
  # w:the length of a sliding window
  # Output: standardized baseline dataframe and residuals dataframe
  dt_ma <- dt %>% group_by(id) %>% 
    mutate(y_ma=TTR::SMA(y,w),res=y-y_ma)
  baseline <- dt_ma %>% na.omit() %>% group_by(id) %>%     
    mutate(y_base=scale(y_ma)) 
  baseline <- as.data.frame(baseline[,colnames(baseline) %in% c('id','days','y_base')])      # apply standardization again
  residual <- dt_ma %>% na.omit() 
  residual <- as.data.frame(residual[,colnames(residual) %in% c('id','days','res')])
  
  return( new("Extract_baseline", 
             base_line = baseline,
             residual = residual))
}
#########

#########
# Caculate SBD distance
SBD_distance <- function(dt,df){
  # shape-based distance on the basis of cross-correlation
  #input:standardized sequence and its length
  ## dt:a standardized sequence
  ## df:a standardized sequence
 
  #output: best slides s,SBD (ranges from 0 to 2)
  n <- length(dt);m <- length(df)  #n:the length of dt, m:the length of df
  RR <- sapply((-n+1):(n-1),function(x){
    if(x >= 0){
      if(m >= n){
        t <- crossprod(dt[(x+1):n],df[1:(n-x)])
      } else {
        t <- crossprod(dt[(x+1):(m+x)],df[1:m])
      }
    } else {
      if (m >= n){
        t <- crossprod(dt[1:(n+x)],df[(m+1-n-x):m])
      } else {
        t <- crossprod(dt[(n-m+x+1):(n+x)],df[1:m])
      }
    }
  })
  NCC <- rbind((-n+1):(n-1),RR/sqrt(crossprod(dt)*crossprod(df)))
  max_s <- NCC[1, which.max(NCC[2,])];NCC_s <- NCC[2, which.max(NCC[2,])]
  SBD <- 1-NCC_s  # ranges from 0 to 2
   return(c(s=max_s,SBD=SBD))
}

# Caculate SBD distance matrix 
# Caculate SBD distance matrix 
SBD_dist <- function(dt,n){
  # Input:dt,n(the number of time series)
  # Output: SBD distance matrix 
  res_mat <- matrix(NA,n,n) 
 
  SBD_mat <- sapply(1:n,function(x){
    i_mat <- dt[,x]
    j_mat <- sapply(x:n,function(y){
      SBD_distance(i_mat,dt[,y])[2]
    })
  })
   for (i in 1:n){
    res_mat[i:n,i] <- SBD_mat[[i]]
  }
  res_mat[upper.tri(res_mat)] <- 0
  r2 <- res_mat+t(res_mat)-diag(diag(res_mat))
  return(r2)
}


# Density-based Clustering:DBSCAN
# Density Estimation
# key parameter:density radius eps
k_distance <- function(mat,minPts,id){
  # Input:
  ## mat:a matrix
  # minPts:a mininum number of the core point's k-Nearest-Neighbor(KNN)
  # Output:
  
  if (anyNA(mat)==TRUE) {
    stop("Please check your matrix with NAs.")
  }
  if (ncol(mat) < 5){
    stop("Data do not have enough points.")
  }
  k_matrix <- as.matrix(mat) 
 
  k_SBD <- apply(k_matrix,2,function(x){
    i_max <- max(sort(x)[2:minPts])
    return(i_max)
  })
  dt_SBD <- data.frame(id=id,SBD=k_SBD)
  
  return(dt_SBD) 
}

# set default
FindCandidateRadius <- function(k_dis,start,end,len_thresh,slope_thresh,slope_diff_thresh) {
  # Input:
  # k_dis:k_dis curve sorted in descending order
  # start:the smallest index that has k_dis[start] <=max_radius
  # end:the last index of k_dis
  # len_thresh:mininum length of a flat part
  # slope_thresh:small value preventing search in steep area
  # slope_diff_thresh:small value indicates a flat part 
  # Output: a list of candidate radiuses
  if ((end-start) <= len_thresh){
    stop("Search area contains few points")
  }
  # set initial number
  r <- -1
  diff <- 2
  for (i in (start+1):(end-1)) {
      leftslope <- (k_dis$SBD[i]-k_dis$SBD[start])/(i-start)
      rightslope <- (k_dis$SBD[end]-k_dis$SBD[i])/(end-i)
      if ((leftslope > slope_thresh) | (rightslope > slope_thresh)){
        print("Search area is steep.")
      } 
      if (abs(leftslope-rightslope) < diff ){
        r <- i
        diff <- abs(leftslope-rightslope)
      }
  }
  
  if (diff < slope_diff_thresh){
    return(r)
  } 
}

######
# dbscan cluster training 
ROCKA_kdist <- function(dt,w=2,minPts=4){
  # Input: 
  # dt :sample data with column id ,days,y 
  ##id ,  days ,      y
  ##1 ,2018-01-01, 100
  # parameters:
  ## minPts:a mininum number of the core point's k-Nearest-Neighbor(KNN)
  #Output:a dataframe with column an object(id) and SBD distance
  
  # Preprocessing
  dt_smooth <- Smooth_anomaly(dt)
  dt_extract <- Extract_Baseline(dt_smooth,w=w) 
  dt_baseline <- dt_extract@base_line  
  
  # Reshape the input dataframe
  baseline_new <- tidyr::spread(dt_baseline,id,y_base) %>% na.omit() # data.frame with column days,id1,id2,...idn
  # Conduct dist object using SBD distance
  base_mat <- as.matrix(baseline_new[,-1])
  nb <- ncol(base_mat)
  #a <- Sys.time()
  tt <- SBD_dist(dt=base_mat,n=nb) # n*n dimension,matrix  baseline_kdis <- k_distance(dt=baseline_new,minPts = minPts)  # takes O(n^2)
  #b <- Sys.time()
  id <- colnames(base_mat)
  baseline_kdis <- k_distance(tt,minPts,id)
  return(baseline_kdis)
}

ROCKA_train <- function(dt,minPts=4,max_radius=0.2,len_thresh=10,slope_thresh=0.01,slope_diff_thresh=0.001){
  # Input: 
  ## dt:a dataframe with column an object(id) and SBD distance
  ## minPt:a mininum number of the core point's k-Nearest-Neighbor(KNN)
  ## max_radius:the largest radius to avoid search in a steep area
  ## len_thresh:mininum length of a flat part
  ## slope_thresh:small value preventing search in steep area
  ## slope_diff_thresh:small value indicates a flat part 
  # Output:gives out an object of class 'dbscan' which is a LIST with components
  
  input_kdis <- dt %>% arrange(desc(SBD)) %>%  
    mutate(idx=seq(1,nrow(dt),1))  # input k_dis to estimate the key parameter:density radius
  
  # compute initial candidated point r
  start_dis <- input_kdis$idx[which.max(input_kdis$SBD < max_radius)]
  end_dis <- nrow(input_kdis)
  r_ini <- FindCandidateRadius(input_kdis,start_dis,end_dis,len_thresh,slope_thresh,slope_diff_thresh) 
  
  # Find all candidates using divide and conquer
  pre_candidates <- data.frame();post_candidates <- data.frame()
  ii <- r_ini;jj <- r_ini
  while(ii > (start_dis+len_thresh+2)) {
      ii=ii-1
      ii <- FindCandidateRadius(input_kdis, start = start_dis, end =ii,len_thresh,slope_thresh,slope_diff_thresh)
      pre_candidates <- rbind(pre_candidates,ii)
    }
  while(jj < (end_dis-len_thresh-2)) {
      jj=jj+1
      jj <- FindCandidateRadius(input_kdis, start = jj, end = end_dis,len_thresh,slope_thresh,slope_diff_thresh)
      post_candidates <- rbind(post_candidates,jj)
    }
   
   r_candidates <- c(as.numeric(pre_candidates[,1]),as.numeric(post_candidates[,1])) 
   r_set <- input_kdis %>% filter(as.numeric(idx) %in% r_candidates) 
   r_best <- r_set$SBD[which.max(r_set$SBD)]
   
   tt1 <- as.dist(tt) # dist object (n-1)*(n-1)
   set.seed(123)
   res <- fpc::dbscan(tt1,eps = r_best, MinPts = minPts,method='dist')
   centroid_dat <- data.frame(id=id,cluster=res$cluster,SBD_sq=rowSums(tt^2)) %>% 
     group_by(cluster) %>% filter(row_number(SBD_sq)==1) %>% as.data.frame()
   return( new("ROCKA_cluster", 
               cluster_res = res$cluster,
               centroid_dat = centroid_dat))
   
  }
######

######
# Assignment:
# assign rest raw time series to the clusters based on the centroids or mark them as outliers
SimilaritytoCentroid <- function(x,y,cluster,w=2){
  # Input:
  ## x: a raw dataframe with column id , days, y 
  ## y: a dataframe of the centroid of each cluster with column id, days, y
  ## w: the length of a sliding window
  # Output:
  ## the SBD distance between raw ts and centroid

  x_smooth <- Smooth_anomaly(x)  # preprocessing
  x_extract <- Extract_Baseline(x_smooth,w=w) # baseline extraction
  x_baseline <- x_extract@base_line  
  
  y_smooth <- Smooth_anomaly(y)  # preprocessing
  y_extract <- Extract_Baseline(y_smooth,w=w) # baseline extraction
  y_baseline <- y_extract@base_line  
  # caculate the SBD between baseline ts and each cluster centroid
  x_in <- data.matrix(x_baseline['y_base'])
  colnames(x_in) <- unique(x_baseline['id'])
  y_in <- data.matrix(y_baseline['y_base'])
  colnames(y_in) <- unique(y_baseline['id'])
  SBD_xy <- SBD_distance(dt=x_in,df=y_in)
  return(c(cluster=as.character(cluster),
           similar_SBD=as.numeric(SBD_xy[2])))
} 

ROCKA_Assignment <- function(r_dt,c_dat,dt,s_thresh=0.2){
  # Input:
  ## r_dt: raw dataframe with column id , days, y
  ## c_dat: centroid dataframe from function ROCKA_train@centroid_dat
  ## dt: raw sample dataset
  ## s_scale: the threshold ,default 0.2.It means a raw ts whose SBD to its nearest cluster centroid 
  ## larger than 0.2 is considered as an outlier.
  ## Output:
  ## the result of assignment
  if(!all(c('id','cluster','SBD_sq') %in% names(c_dat))) {
    stop("Data have no columns 'id' or 'cluster' or 'SBD_sq'. Please check your data.")
  }
  c_dat$id <- as.character(c_dat$id)
  n <- nrow(c_dat)
  s <- 2
  cluster_i <- NULL
  for (i in 1:n){
    yi <- dt %>% filter(id==as.numeric(c_dat$id[1]))
    si <- SimilaritytoCentroid(x=r_dt,y=yi,cluster=c_dat$cluster[i],w=2)
    if (as.numeric(si[2]) <= s){
      s <- as.numeric(si[2])
      cluster_i <- c_dat$cluster[i]
    }
  }
  # raw time serie whose SBD to its nearest cluster centroid larger than 0.2 is considered as an outlier
  if (s <= as.numeric(s_thresh)){
    assign_c <- c(id=as.character(r_dt$id[1]),cluster=as.character(cluster_i),similarity=s)
  } else {
    assign_c <- c(id=as.character(r_dt$id[1]),cluster="outlier",similarity=s)
  }
  
  return(assign_c)
}


######






















