#The goal of the script is to automate the analysis of plant vascular system in cross-section

#More than simply providing metrics related to xylem vessel diameter, it is also providing metrics
#related to vessels spatial arrangement and connectivity while considering their true shapes in cross-section

#The script can handle large samples such as whole stem cross-section having thousands of vessels in few seconds
#Well, the kind of samples you couldn't handle using regular procedure   

#The following metrics are return:
#Dh = Hydraulic diameter
#Kh_th = Theoretical hydraulic conductivity (kg.m.s-1.MPa-1)
#ks_th = Theoretical sapwood specific hydraulic conductivity (kg.m-1.s-1.MPa-1)
#F = Lumen fraction
#Fc = Contact fraction
#Vessel Grouping Indices = GI 
#Vd = Vessel density

#and more could be done for specific needs but it's time to share

#Please see Scholz, A., Klepsch, M., Karimi, Z. and Jansen, S., 2013. How to quantify conduits in wood?. Frontiers in plant science, 4, p.56.

#The input is a .txt file one can obtain with ImageJ (FIJI is preferred) from the "Bare Outlines" picture produced
#by the "Analyse Particules" tool
#X and Y coordinates are then exported using the  Analyze > Tools > Save XY Coordinates...
#to create the .txt file

#Loading dbscan package
library(dbscan)

#enter the name of ".txt" file here
#Such file is obtained with ImageJ from the "Bare Outlines" picture produced by the "Analyse Particules" tool
#then export X and Y coordinates using the  Analyze > Tools > Save XY Coordinates... in ImageJ (FIJI is preferred) to create the .txt file
Sample = "Drawing"
#creating the path to load the data
#here the demo_set.txt is just on "D/:" but it might differ for you
path <- paste("D:/",paste(Sample,".txt", sep = ""), sep = "")
#creating a dataframe where our data will be loaded
Dataset <- read.table(path,header=TRUE, stringsAsFactors=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) 

#let's plot X and Y coordinates and see what it looks like
#bare outlines of xylem vessels...  
plot(Dataset[,c(1,2)], pch=1, cex=0.01)

#now it is needed to enter some parameters
#The one given hereafter were determined for grapevine in my experimental setting 

#Xa is the xylem area analyzed in mm2 (Ax = 1.236521 is my demo_set)
Ax = 1.236521
#CWmin is the theoretical minimum distance between vessels (in pixel) (CWmin = 3 in my case)
#CWmin is the intercept of the regression line between the double vessel wall thichness (Tvw)
#and the sum of the diameter of vessel being in contact      
CWmin = 3
#TvwI is the slope of the regression described above, 0.05 in the plant species (Vitis vinifera) displayed here
#These parameters will be used further on to weight the distance
#at which a contact between vessels will be detected
TvwI=0.05
#the Scale is the number of pixels per microns, and depends on your experimental settings, mine is 0.5507
Scale=0.5507
#Cl1 is the distance (in pixels) we will use to extract clusters of pixels
#corresponding to individual vessels from the dbscan analysis
#Cl1=3 should work in any case 
Cl1 = 3
#Cl2 is the distance (in pixels) we will use to extract clusters of pixels
#corresponding to groups of vessels with dbscan
#Cl2 might depend on the species, picture resolution, and so on...
#This should be a bit above the largest double vessel wall thichness (Tvw) measured in your samples     
Cl2 = 15

#Creating the dataframe with X and Y coordinates of pixels
map <- Dataset[ ,c(1,2)] 
#Creating the clusters with the dbscan algorithm
res <- optics(map, eps = 20,  minPts = 5) 
#Extracting clusters of pixels from res with the distance Cl1 = individual vessels
res1 <- extractDBSCAN(res, eps_cl =Cl1) 
#Extracting cluster of pixels from res with the distance of Cl2 = groups of vessels
res2 <- extractDBSCAN(res, eps_cl = Cl2)

#Creating an empty vector where the vessel area to which the pixels belong (V_area) will be stored
V_area<-vector("numeric",NROW(map)) 
#Creating an empty vector where the vessel diameter to which the pixels belong (V_diameter) will be stored
V_diameter<-vector("numeric",NROW(map)) 
#Creating an empty vector where pixels in contact will be labelled (Pix_contact)
Pix_contact<-vector("numeric",NROW(map)) 
#Creating a dataframe (df1) with map, cluster’s information, internal sequence of pixels, and the empty columns for V_area, V_diameter and Pix_contact to be stored
df1<-cbind(map, res1$cluster, res2$cluster,seq.int(NROW(map)),V_area,V_diameter, Pix_contact) 


#Let's check what our groups of vessels looks like
df1_fact<-df1
colnames(df1_fact)[4]<-"VG"
as.factor(df1_fact[,4])
df1_agg<-aggregate(df1_fact, list(df1_fact[,4]), FUN=mean, simplify = TRUE,drop = TRUE, discrete=F)
plot(df1[,c(1,2)], col=df1[,4], pch=1, cex=0.01)
text(df1_agg[,2],df1_agg[,3],df1_agg[,1])
#Victory!?
#Not really... look at cluster 27, 21, 12, 24. dbscan detected larger vessel's groups then it should
#But this is what we want here, because it can be fixed later on within each group
#using another strategy (hclust) while weighting the Tvw we expect between vessel's pairs, based on their diameters


#The nested loop below aim at computing the area of vessels
#It takes each vessel one by one based on the values from res1
#It sums the number of pixels in each line of pixels the vessel is composed = vessel’s area
#It returns the values to where it belongs for each pixel in df1$V_area based on the values of res1
for (i in min(df1[,3]):max(df1[,3]))
{
  V<-df1[which(df1[,3]==i),]
  area<-vector()
  for (j in min(df1[which(df1[,3]==i),1]):max(df1[which(df1[,3]==i),1]))
  {
    point_line<-V[which(V[,1]==j),]
    line<-max(point_line[,2])-min(point_line[,2])+1
    line<-ifelse(line<=0,0,line)
    area[j]<-line
  }
  Varea<-sum(area, na.rm = TRUE)
  df1$V_area[which(df1[,3]==i)]<-Varea
}

res3<-numeric(length=NROW(df1))
res3<-as.data.frame(res3)
res3<-cbind(res3, seq.int(NROW(res3)))

#This function is meant to create matrix of the sum of a given parameters x
summat = function(x){
  D = matrix(as.numeric(NA), NROW(x), NROW(x))
  for (j in 1:NROW(x)){
    d = x[[j]] + x[-j]
    D[j,-j] = d
  }
  if (!all(is.na(diag(D)))){
    stop("Not all diagonal elements zero")
  }
  diag(D) = 0
  if (!is.null(names(x))) colnames(D) = rownames(D) = names(x)
  return(D)
}

i=1
#This loop will deal with vessel's groups (found with dbscan) one by one
#It does detect pixels below a weighted distance with the vessel's group and return it to df1$Pix_contact column
#Then it performs an hclust on a matrix which is the difference between the distance matrix of pixels and
#the expected Tvw of vessel's pairs (weighted by the diameter of vessels)
#The hclust tree is cut at distance = 0, meaning the clusters returned (i.e. vessel's groups)
#are for pixel distance values lower than the expected Tvw   
nr=NROW(table(df1[,4]))
for (i in 1: nr)
{
  #Creating a small dataset with the data from the dbscan vessel's groups   
  map2<-df1[df1[,4]==i,]   
  #Creating a distance matrix with the coordinates of pixels in map2
  dist_map <- dist(map2[,c(1,2)], method = "euclidean", diag = FALSE, upper = TRUE)
  #Converting the distance matrix into matrix
  dist_map <- as.matrix(dist_map)
  #Computing the theoretical diameter of vessels
  D_vessel<-2*sqrt(map2[,6]/pi)
  #Computing the expected weighted thicknesses of the vessel's walls
  CW_pred<-TvwI*D_vessel+CWmin
  #Creating a matrix of the expected double vessel wall thickness of vessel's pairs (Tvw)
  CW_pred_Matrix<-summat(CW_pred)
  #Creating a logical distance matrix where pixels from the same vessel = 0 
  logic<-dist(map2[,3], method = "euclidean", diag = FALSE, upper = TRUE)
  #Converting the distance matrix into matrix
  logic<-as.matrix(logic)
  #Turning to 0 the pixel's distances that are larger than the expected Tvw 
  d2<-ifelse(dist_map<=CW_pred_Matrix,dist_map,0)
  #Turning to 0 distance when they were computed within the same vessel 
  d3<-d2*logic
  #Returning the max value per line in the d3 matrix 
  kept<-apply(d3, 1, max, na.rm=TRUE)
  #Converting kept in vector
  kept<-as.vector(kept)
  #Converting all values that are not equal to 0 to 1, 1 meaning pixels belong to the intervessel walls
  kept<-ifelse((kept)!=0,1,0)
  #Adding the kept vector to map2
  map2<-cbind(map2, kept)
  #Storing the identified pixels in res3
  res3[intersect(map2[,5],res3[,2]),1] <- map2[,9]
  
  #Here is the hclust procedure to re-analyze the vessel's groups with the expected weighted Tvw  
  #Substracting the matrix of expected Tvw from the distance matrix of vessels 
  Contact_matrix <- dist_map-CW_pred_Matrix
  #Converting the matrix into distance matrix object so that hclust can run
  Contact_matrix <- as.dist(Contact_matrix)
  #Running hclust in single mode
  hc <- hclust(Contact_matrix, method= "single")
  #Cutting the hclust tree at distance = 0 
  hc_cut <- cutree(hc, h=0)
  #Adding the hc_cut vector to map2
  map2<-cbind(map2, hc_cut)
  #Storing the identified clusters in res3
  res3[intersect(map2[,5],res3[,2]),3] <- map2[,10]
}

#Returning the results of Pix_contact in df1 
df1$Pix_contact <- res3[,1]
#Returning the results of the hclust in df1
df1$hc_cut <- res3[,3]
#Creating a "compound" variable for vessel's groups with the new (hclust) and old (dbscan) analysis outputs
#if an old (dbscan) vessel's group x was splitted in 2 or n vessel's groups so now they are named x.1, x.2, ..., x.n
#This should be problematic if n>9, chech Cl2 parameters if needed. 
df1$cons_clust <- df1$hc_cut*0.1+df1[,4]

#let's check what our groups of vessels looks like now!
df1_fact<-df1
colnames(df1_fact)[10]<-"VG"
df1_agg<-aggregate(df1_fact, list(df1_fact$VG), FUN=mean, simplify = TRUE,drop = TRUE)
plot(df1[,c(1,2)], col=df1[,4], pch=1, cex=0.01)
text(df1_agg[,2],df1_agg[,3],df1_agg[,1])
#Looks pretty accurate to me!


#let's check our vessel contact detection
#Neat!
plot(map[,c(1,2)], col=df1$Pix_contact+1, pch=1, cex=0.1+0.5*df1$Pix_contact)

#Now we can compute all metrics
#Calculating the number of vessels from the count of different number in res1
Nvessel<-NROW(table(df1[,3]))
#Calculating the number of vessels from the count of different number in res1
Ngroup<-NROW(table(df1$cons_clust))
#Calculating the grouping indices of the sample 
GI<-Nvessel/Ngroup
#Calculating the diameter of xylem vessels
df1$V_diameter<-(2*sqrt(df1$V_area/pi))/Scale
#Calculating the hydraulic diameter of the sample
Dh<-sum((tapply(df1$V_diameter,df1[,3], mean))^5)/sum((tapply(df1$V_diameter,df1[,3], mean))^4)
#Calculating the theoretical hydraulic conductivity Kh_th (kg.m.s-1.MPa-1)
Kh_th<-1000*sum(pi*((10^(-6))*(tapply(df1$V_diameter,df1[,3], mean)))^4)/(128*(1.002*10^(-9)))
#Calculating the theoretical sapwood specific hydraulic conductivity Ks_th (kg.m-1.s-1.MPa-1)
Ks_th<-10^6*(Kh_th/Ax)
#Calculating the sum of the area of vessels (in mm2)
sum_Varea<-(sum((pi*(tapply(df1$V_diameter,df1[,3], mean)/2)^2)))/10^6
#Calculating the vessel lumen fraction F
F <- sum_Varea/Ax 
#Calculating the contact fraction Fc
Fc <- sum(df1[,8])/NROW(df1)
#Calculating the vessel density Vd (count / mm2)
Vd <- Nvessel/Ax

#Exporting the analysis output of vessel's groups as a .tiff figure
df1_fact<-df1
colnames(df1_fact)[10]<-"VG"
df1_agg<-aggregate(df1_fact, list(df1_fact$VG), FUN=mean, simplify = TRUE,drop = TRUE)
tiff(paste(Sample,"_groups.tiff", sep = ""), units="in", width=10, height=10, res=600, compression = 'lzw')
plot(df1[,c(1,2)], col=df1[,4], pch=1, cex=0.01)
text(df1_agg[,2],df1_agg[,3],df1_agg[,1])
dev.off()

#Exporting the analysis output of intervessel wall as a .tiff figure
tiff(paste(Sample,"_contact.tiff", sep = ""), units="in", width=10, height=10, res=600, compression = 'lzw')
plot(map[,c(1,2)], col=df1$Pix_contact+1, pch=1, cex=0.1+0.5*df1$Pix_contact)
dev.off()

#Creating an data.frame were metrics are stored
output<-cbind(Dh, Kh_th, Ks_th, Ax, Vd, sum_Varea, F, Fc, Nvessel, Ngroup, GI)
rownames(output) <- Sample
output

#Enjoy!

#Jerome Pouzoulet, PhD