## data preparation

library(xlsx)
library(ggmap)
library(mapproj)

setwd("~/Desktop/ploy_match")
data<-read.xlsx(file="sample_delarea_coordinates2.xlsx",1,header=T)

test <- data[,2]
loca_pattern <- "[0-9]+.[0-9]+, [0-9]+.[0-9]+"
matches <- gregexpr(pattern = loca_pattern, text = test)
loca<-regmatches(test, matches)

matlist <- NULL
for (i in 1:length(loca))
{
  n <- length(loca[[i]])
  matlist[[i]]<-matrix(as.numeric(unlist(strsplit(loca[[i]],", "))),n,2,byrow = T)
  matlist[[i]]<-matlist[[i]][-1,]
}

## the 1st,5th and 7th are reversed
matlist[[1]]<-cbind(rev(matlist[[1]][,1]),rev(matlist[[1]][,2]))
matlist[[5]]<-cbind(rev(matlist[[5]][,1]),rev(matlist[[5]][,2]))
matlist[[7]]<-cbind(rev(matlist[[7]][,1]),rev(matlist[[7]][,2]))


## function


calculate_distance<-function(mat1,mat2,start)
{
  
  ## rearrange the mat
  
  
  if (start != 1)
    mat1 <- rbind(mat1[start:nrow(mat1),],mat1[1:(start-1),])
  
  mat1 <- rbind(mat1,mat1[1,])
  mat2 <- rbind(mat2,mat2[1,])
  
  list <- list(mat1,mat2)
  
  
  ## calculate theta-function
  
  for (j in 1:2)
  {
    
    ## calculate the angles,angle rotates,concavity
    ## prepare for calculating theta-function
    
    mat <- list[[j]]
    k <- rep(NA,nrow(mat))
    theta <- rep(NA,nrow(mat))
    crossp <- rep(NA,nrow(mat)-1)
    dist <- rep(NA,nrow(mat)-1)
    for (i in 1:(nrow(mat)))
    {
      if (i>1 & i<nrow(mat))  
      {
        x<- c((mat[i,1] - mat[i-1,1]) , (mat[i,2] - mat[i-1,2]))
        y<- c((mat[i+1,1] - mat[i,1]) , (mat[i+1,2] - mat[i,2]))
        dist[i]<-dist(rbind(mat[i,],mat[i+1,]),method = "euclidean")
        crossp[i]<- x[1]*y[2]-x[2]*y[1]
      }
      else if (i==nrow(mat))  
      {
        x<- c((mat[i,1] - mat[i-1,1]) , (mat[i,2] - mat[i-1,2]))
        y<- c((mat[2,1] - mat[1,1]) , (mat[2,2] - mat[1,2]))
      }
      else
      {
        x<- c(1,0)
        y<- c((mat[i+1,1] - mat[i,1]) , (mat[i+1,2] - mat[i,2]))
        x1<- c((mat[nrow(mat),1] - mat[nrow(mat)-1,1]) , (mat[nrow(mat),2] - mat[nrow(mat)-1,2]))
        y1<- c((mat[2,1] - mat[1,1]) , (mat[2,2] - mat[1,2]))
        dist[i]<-dist(rbind(mat[i,],mat[i+1,]),method = "euclidean")
        crossp[i]<- x1[1]*y1[2]-x1[2]*y1[1]
      } 
      k[i] <- acos(sum(x * y) / sqrt(sum(x^2) * sum(y^2)))
      
      ## calculate theta-function  
      
      rotate <- ifelse(crossp>0, -k ,k)
      theta[1]<-0
      for(i in 2:nrow(mat))
        theta[i] <- rotate[i-1]+theta[i-1]
      theta <- theta[-1]
    }
    if(j==1)
      theta_s1 <- rbind(dist,theta)
    else
      theta_s2 <- rbind(dist,theta)
  }
  
  #print(theta_s1)
  #print(theta_s2)
  
  ## rescale
  theta_s2[1,]<-theta_s2[1,] * sum(theta_s1[1,])/sum(theta_s2[1,])
  
  ## cutting into rect strips
  
  list<-list(theta_s1[1,],theta_s2[1,])
  list2<-list(theta_s1[2,],theta_s2[2,])
  
  new_s<-NULL
  theta1<-NULL
  theta2<-NULL
  
  i=1 ## for list[[1]]
  j=1 ## for list[[2]]
  k=1 ## for new_s
  theta1[k]<-list2[[1]][i]
  theta2[k]<-list2[[2]][j]
  new_s[k]<-min(list[[1]][i],list[[2]][j]) ## choose the minimal one to cut
  left_index<-3-which.min(c(list[[1]][i],list[[2]][j]))  ## choose the rest one that has distance left
  if(left_index == 1)
  {
    left<-list[[left_index]][i] - new_s[k]
    j=j+1
  }else
  {
    left<-list[[left_index]][j] - new_s[k]
    i=i+1 
  }
  k=k+1
  winner <- 3-left_index
  times <- 1
  
  while( ifelse(left_index == 1,length(list[[3-left_index]]) >= j,length(list[[3-left_index]]) >= i))
  {
    theta1[k]<-list2[[1]][i]
    theta2[k]<-list2[[2]][j]
    if(left_index ==1)
    {
      new_s[k]<-min(list[[3-left_index]][j],left) ## choose the minimal one to cut
      ## 3-left_index is the winner of last turn
      ## the following left_index is the loser of this turn
      left_index<-ifelse(which.min(c(list[[3-left_index]][j],left))==1,left_index,3-left_index) 
    }
    else 
    {
      new_s[k]<-min(list[[3-left_index]][i],left)
      left_index<-ifelse(which.min(c(list[[3-left_index]][i],left))==1,left_index,3-left_index) 
    }
    
    if (winner == 3-left_index) 
    {times = times + 1} else
    {
      winner = 3-left_index
      times = 1
    }
    
    if(left_index == 1)
    {
      left<-list[[left_index]][i] - sum(new_s[(k-times+1):k])
      j=j+1
    }
    else
    {
      left<-list[[left_index]][j] - sum(new_s[(k-times+1):k])
      i=i+1 
    }
    k=k+1
  }
  y <- theta1^2 - 2*theta1*theta2 + theta2^2
  return(sum(new_s * y))
  
}


## use the function to get the result
mat <- matrix(NA,7,7)
for (i1 in 1:7)
{
  for (i2 in i1:7)
  {
    result <- NULL
    for (i in 1:nrow(matlist[[i1]]))
    { result[i]<-calculate_distance(matlist[[i1]],matlist[[i2]],start=i) }
    
    ## final_result
    mat[i1,i2]<-10*min(result)
    #print(paste(i1,"vs",i2,":",min(result)))
  }
}

write.csv(mat,file="result.csv")


## if you want to draw the ith map shape
## take the first one as an example

i<-1
data <- matlist[[i]]
colnames(data)<-c("long","lat")
lonc <- apply(data,2,mean)[1]
latc <- apply(data,2,mean)[2]
ggmap(get_googlemap(center=c(lonc,latc), zoom=13, maptype='roadmap'), extent='device') + geom_polygon(data = as.data.frame(data), aes(x=long, y = lat), fill = NA, color = "red") + coord_fixed(1.3) 

