library(raster)
library(Matrix)
library(dplyr)
library(parallel)
options(future.globals.maxSize = 4000 * 1024^2)

labelling<- function(mat){
  tmp<-data.frame(matrix(ncol=3, nrow=0))
  colnames(tmp)<-c("x","y","name")
  for(i in (1:max(mat))){
    ind<-center_spark(mat,i)
    tmp<-rbind(tmp, c(ind[2],nrow(mat)-ind[1], i))
  }
  return(tmp)
}
affichage <- function(img_mat_fun,mat){
  res_img<-matrix(0,dim(img_mat_fun)[1],dim(img_mat_fun)[2])
  if(max(mat)!=0){
    for(i in (1:max(mat))){
      ind<-which(mat==i, arr.ind=TRUE)
      for(el in 1:nrow(ind)){
        res_img[ind[el,1],ind[el,2]]<-img_mat_fun[ind[el,1],ind[el,2]]
      }
      
    }
  }
  return(res_img)
}


center_spark<-function(mat_spark, sparknumber){
  ind<-which(mat_spark==sparknumber, arr.ind=TRUE)
  x<-mean(ind[,1])
  y<-mean(ind[,2])
  return(c(x,y))
}
results_dir<-paste0("results/", file_name)
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}


img=raster(paste0("data/",file_name,".tif"))
img_mat=as.matrix(img)

filtre<- function(pixel){
  if (pixel<(moyenne+85)){
    return(0L)
  }else{
    return(1L)
  }
}

moyenne=mean(img_mat)
img_mat_ft=matrix(unlist(mclapply(img_mat, filtre, mc.cores=3)), nrow=dim(img_mat)[1], ncol=dim(img_mat)[2])
img_mat_ft[dim(img_mat_ft)[1],dim(img_mat_ft)[2]]=1

find_transient<- function(mat){
  transient_x_coord<-c()
  for(ix in 1:nrow(mat)){
    if(sum(mat[ix,]!=0)/ncol(mat)>0.5){
      transient_x_coord<-c(transient_x_coord,ix)
    }
  }
  transients<-c()
  for(el in transient_x_coord){
    transients<-c(transients, el:(min(el+30,nrow(mat))))
  }
  transients<-unique(transients)
  return(transients)
}



find_blanks<- function(mat){
  blanks<-c()
  for(ix in 1:nrow(mat)){
    if(sum(mat[ix,]!=0)/ncol(mat)<0.16){
      blanks<-c(blanks,ix)
    }
  }
  return(blanks)
}



transient <- find_transient(img_mat_ft)


img_mat_ft[transient,]=0
blanks<-find_blanks(img_mat_ft)
blanks<-sort(unique(c(1,nrow(img_mat_ft),blanks)))

while(max(blanks[2:length(blanks)]-blanks[1:length(blanks)-1])>200){
  ind<-which.max(blanks[2:length(blanks)]-blanks[1:length(blanks)-1])
  ind2<-which.min(rowSums(img_mat_ft[(blanks[ind]+1):(blanks[ind+1]-1),]))
  blanks<<-c(blanks,(ind2+blanks[ind]))
  blanks<<-sort(blanks)
  print((ind2+blanks[ind]))
}


analyse<-function(img_mat_ft_fun,img_mat_fun,renum2){
  img_mat_ft_fun[dim(img_mat_ft_fun)[1],dim(img_mat_ft_fun)[2]]=1
  img_mat_ft_fun <- Matrix(img_mat_ft_fun, sparse = TRUE)  
  
  sm = summary(img_mat_ft_fun)
  d <- dist(sm, "manhattan")
  if(!(dim(as.matrix(d))[[1]]<=1 | dim(as.matrix(d))[[2]]<=1)){
    gr = cutree(hclust(d, "single"), h = 1)
    mat_sparks=sparseMatrix(i = sm[, "i"], j = sm[, "j"], x = gr)
    rm(gr)
  }else {
    mat_sparks<- Matrix(0,nrow=nrow(img_mat_ft_fun),ncol=ncol(img_mat_ft_fun), sparse = TRUE)
  }

  rm(sm)
  rm(d)
  
  renum=1
  for(i in (1:max(mat_sparks))){
    if(sum(mat_sparks==i)<40 ){
      mat_sparks[mat_sparks==i]=0
    }else{
      mat_sparks[mat_sparks==i]=renum
      renum=renum+1
    }
  }
  mat_sparks=as.matrix(mat_sparks)
  
  
  
  
  
  
  
  
  distance_pixel<-function(x1,y1,x2,y2){
    distance<-sqrt((x1-x2)^2+(y1-y2)^2)
    return(distance)
  }
  
  
  subspark<- function(mat, sparknumber, mat_img_fun){
    mat[mat!=sparknumber]=0
    mat[mat==sparknumber]=1
    spark<-affichage(mat_img_fun, mat)
    moy_sp=mean(spark)
    std_sp=sd(spark)
    filtre_sp<- function(pixel){
      if (pixel>(moy_sp+3*std_sp)){
        return(1L)
      }else{
        return(0L)
      }
    }
    peaks=matrix(unlist(mclapply(spark, filtre_sp, mc.cores=3)), nrow=dim(spark)[1], ncol=dim(spark)[2])
    peaks=Matrix(peaks, sparse = TRUE)
    sm = summary(peaks)
    d <<- dist(sm, "manhattan")
    if(!(dim(as.matrix(d))[[1]]<=1 | dim(as.matrix(d))[[2]]<=1)){
      gr = cutree(hclust(d, "single"), h = 1)
      res=sparseMatrix(i = sm[, "i"], j = sm[, "j"], x = gr)
      rm(gr)
    }else {
      res<- Matrix(0,nrow=nrow(peaks),ncol=ncol(peaks), sparse = TRUE)
    }

    rm(d)
    rm(sm)
    renum=1
    for(i in (1:max(res))){
      if(sum(res==i)<15 ){
        res[res==i]=0
      }else{
        res[res==i]=renum
        renum=renum+1
      }
    }
    res=as.matrix(res)
    if(max(res)==0){
      return(mat)
    }else{
      newspark_center<-list()
      for(newspark in 1:max(res)){
        newspark_center[length(newspark_center)+1]<-list(center_spark(res,newspark))
      }
      inds<-which(mat==1, arr.ind=TRUE)
      for(i in 1:nrow(inds)){
        distances<-c()
        for(newspark in 1:max(res)){
          distances<-c(distances,(distance_pixel(newspark_center[[newspark]][1],newspark_center[[newspark]][2],inds[i,1],inds[i,2])))
        }
        mat[inds[i,1],inds[i,2]]<-which(distances==min(distances), arr.ind=TRUE)
      }
      return(mat)
    }
  }
  
  
  mat_subsparks<-matrix(0,dim(mat_sparks)[1],dim(mat_sparks)[2])
  if(max(mat_sparks)!=0){
    for(spark in 1:max(mat_sparks)){
      ressub<-subspark(mat_sparks,spark,img_mat_fun)
      for(i in 1:max(ressub)){
        ind<-which(ressub==i, arr.ind=TRUE)
        renum2=renum2+1
        for(el in 1:nrow(ind)){
          mat_subsparks[ind[el,1],ind[el,2]]<-renum2
          
        }
      }
    }
  }
  
  return(list(mat_subsparks,renum2))
}

res_assemblee<-matrix(0,dim(img_mat_ft)[1],dim(img_mat_ft)[2])
if(blanks[1]!=1){
  temp<-analyse(img_mat_ft[1:blanks[1],],img_mat[1:blanks[1],],0)
  res_assemblee[1:blanks[1],]<-temp[[1]]
  
}else{
  temp<-c(0,0,0)
}

for(i in 1:(length(blanks)-1)){
  print(i)
  print((length(blanks)-1))
  if((blanks[i]+1)!=(blanks[i+1])){
    temp<-analyse(img_mat_ft[blanks[i]:blanks[i+1],],img_mat[blanks[i]:blanks[i+1],],temp[[2]])
    res_assemblee[blanks[i]:blanks[i+1],]<-temp[[1]]
  }
}

renum=1
for(i in (1:max(res_assemblee))){
  if(sum(res_assemblee==i)>1000 ){
    res_assemblee[res_assemblee==i]=0
  }else{
    res_assemblee[res_assemblee==i]=renum
    renum=renum+1
  }
}


# resa = raster(res_assemblee)
# extent(resa) = c(1,ncol(res_assemblee),1,nrow(res_assemblee))
# plot(resa)
# textimg_mat<-labelling(res_assemblee)
# colnames(textimg_mat)<-c("x","y","name")
# text(x = textimg_mat$x, y = textimg_mat$y, labels = textimg_mat$name,cex=0.5)





abc<-affichage(img_mat,res_assemblee)
abc[transient,]=0
if(max(res_assemblee)!=0){
  try({
    png(paste0(results_dir,"/plot.png"),height=nrow(abc)*0.3,width = ncol(abc), unit="px")
    resa=raster(abc)
    extent(resa) = c(1,ncol(abc),1,nrow(abc))
    par(mar = c(1, 1, 1, 1)) 
    plot(resa)
    textimg_mat<-labelling(res_assemblee)
    colnames(textimg_mat)<-c("x","y","name")
    text(x = textimg_mat$x, y = textimg_mat$y, labels = textimg_mat$name)
    
    dev.off()
  })
}


writeRaster(resa,paste0(results_dir,"/res_assemblee.tif"), overwrite=TRUE)
sparks_params<-as.data.frame(matrix(ncol = 15))
colnames(sparks_params)<-c("taille","duree", 'debut','fin', 'distance','Transitoire avant','Transitoire après','distance_entre_transitoire','Distance normalisée','surface','f','f0','fsurf0', 'fsurf0v2','Profile Plot')
if(max(res_assemblee)!=0){
  for(i in 1:max(res_assemblee)){
    ind<-which(res_assemblee==i, arr.ind=TRUE)
    debut<-min(ind[,1])
    fin<-max(ind[,1])
    duree<-(max(ind[,1])-min(ind[,1]))*1.88
    taille<-(max(ind[,2])-min(ind[,2]))*0.26
    x<-center_spark(res_assemblee,i)[[1]]
    tmp<-transient-x
    distance<-min(tmp[tmp>=0])*1.88
    transient_apres<-min(tmp[tmp>=0])+x
    transient_avant<-max(tmp[tmp<=0])+x
    distance_entre_transitoire<-(transient_apres-transient_avant)*1.88
    distance_norm<-distance/(transient_apres-transient_avant)
    surface<-length(which(res_assemblee==i))
    profile_plot<-c()
    all_values_of_spark<-c()
    for(i2 in 0:(max(ind[,1])-min(ind[,1]))){
      ind2<-filter(as.data.frame(ind), row==(min(ind[,1])+i2))
      tmp2<-0L
      for(i3 in ind2[,2]){
        tmp2<-tmp2+abc[(min(ind[,1])+i2),i3]
        all_values_of_spark<-c(all_values_of_spark,abc[(min(ind[,1])+i2),i3])
      }
      profile_plot<-c(profile_plot,tmp2/length(ind2[,2]))
    }
    f<-max(profile_plot[!is.nan(profile_plot)])
    f0<-min(profile_plot[!is.nan(profile_plot)])
    fsurf0<-f/f0
    fsurf0v2<-mean(sort(all_values_of_spark, decreasing = TRUE)[1:5])/mean(sort(all_values_of_spark)[1:5])
    sparks_params[i,]<-c(taille,duree,debut,fin,distance,transient_avant,transient_apres,distance_entre_transitoire,distance_norm,surface,f,f0,fsurf0,fsurf0v2, 0)
    sparks_params$'Profile Plot'[i]<-list(profile_plot)
  }
}


# plot(sparks_params$taille)
# plot(sparks_params$duree)
# plot(sparks_params$distance)
# plot(sparks_params$surface)
# plot(1:length(sparks_params$'Profile Plot'[[1]])*1.88,sparks_params$'Profile Plot'[[1]], type='b')

saveRDS(sparks_params, paste(results_dir,"/sparks_params.rds",sep=""))
saveRDS(res_assemblee, paste(results_dir,"/res_assemblee.rds", sep=""))
write.csv(x = subset(sparks_params, select=-`Profile Plot`), file = paste(results_dir, "/sparks_params.csv", sep=""))
conn<-file(paste(results_dir,"/resultats.txt", sep=""))
data_to_write<-c(
  paste("Nombre de sparks", max(res_assemblee)),
  paste("Nb sparks normalisé pour 1s et 20µm", max(res_assemblee)*1000*1.88/dim(res_assemblee)[1]*20*0.26/dim(res_assemblee)[2]),
  paste("Moyenne de la taille (µm)", mean(sparks_params$taille[!is.infinite(sparks_params$taille)])),
  paste("Moyenne de la durée (ms)", mean(sparks_params$duree[!is.infinite(sparks_params$duree)])),
  paste("Moyenne de la distance aux transients (si présent)(ms)", mean(sparks_params$distance[!is.infinite(sparks_params$distance)])),
  paste("F/F0 ", mean(sparks_params$fsurf0)),
  paste("F/F0v2 ", mean(sparks_params$fsurf0v2)),
  paste("Nb sparks dist <20 ", nrow(sparks_params %>% filter(distance < 20))),
  paste("Nb sparks dist <30 ", nrow(sparks_params %>% filter(distance < 30))),
  paste("Nb sparks dist <40 ", nrow(sparks_params %>% filter(distance < 40))),
  paste("Nb sparks dist/dist entre transitoire <0.1 ", nrow(sparks_params %>% filter((distance/distance_entre_transitoire) < 0.1))),
  paste("Nb sparks dist/dist entre transitoire <0.2 ", nrow(sparks_params %>% filter((distance/distance_entre_transitoire) < 0.2))),
  paste("Nb sparks dist/dist entre transitoire <0.3 ", nrow(sparks_params %>% filter((distance/distance_entre_transitoire) < 0.3)))
)          
writeLines(data_to_write, conn)
close(conn)
