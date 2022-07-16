GI1<-function(X,k){          # Gabaix & Ibragimov - reversed A&M
  X<-sort(X)
  n<-length(X)
  if(missing(k)) k<-round(sqrt(n))
  
  x<--log((c(n:1)-0.5)/n)[(n-k+1):n]
  y<-log(X[(n-k+1):n])
  f<-lm(x~y)
  f$coefficients[2]
}

Hill<-function(X,k){          #Hill estimator
  X<-sort(X)
  n<-length(X)
  if(missing(k)) k<-round(sqrt(n))
  (sum(log(X[(n-k+1):n]))/k-log(X[n-k]))^(-1)
}


library(ggplot2)
library(ggeasy)

#names(LogLogWork)
#varlist <- names(LogLogWork)[2:25]
library(readxl)
dataa <- read.csv("covid_rate_confirmed_cases_columns_reordered_March_2021.csv")
varlist <- c('United Kingdom',
             'Germany',
             'France',
             'Italy',
             'Spain',
             'Russia',
             'Netherlands',
             'Sweden',
             'India',
             'Austria',
             'Finland',
             'Ireland',
             'US',
             'US',
             'Lithuania',
             'Canada',
             'Brazil',
             'Mexico',
             'Argentina',
             'Japan',
             'China',
             'Korea',
             'Indonesia',
             'Australia',
             'Australia')

for(i in 1:24){
  y <- dataa[,i*4-1]
  y <- y[y>0]
  x <- y[2:length(y)]-y[1:(length(y)-1)] # first difference
  z <- x[2:length(x)]-x[1:(length(x)-1)] # second difference
  




#for (i in 2:(length(varlist)+1)){
  
N<-length(z)#sum(LogLogWork[!is.na(LogLogWork[,i]),i]>0)
N
Matr2<-replicate(floor(0.15*N)-floor(0.025*N)+1, 0)
SE<-replicate(floor(0.15*N)-floor(0.025*N)+1, 0)
CILeft<-replicate(floor(0.15*N)-floor(0.025*N)+1, 0)
CIRight<-replicate(floor(0.15*N)-floor(0.025*N)+1, 0)
for (k in floor(0.025*N):floor(0.15*N)){
  #lin1<-lm(LogLogWork$lRank[1:k]~as.matrix(LogLogWork[1:k,i]), data=LogLogWork)
  C<-Hill(z,k)#GI1(z,k)#coef(lin1)
  Matr2[k-floor(0.025*N)+1]<- C
  SE[k-floor(0.025*N)+1]<- C*sqrt(2)/sqrt(k)
  CILeft[k-floor(0.025*N)+1]<- C-1.96*C/sqrt(k)
  CIRight[k-floor(0.025*N)+1]<- C+1.96*C/sqrt(k)
}

k1<-floor(0.025*N):floor(0.15*N)
k1 <- k1/N*100
Matr2  
SE
CILeft
CIRight
df <- data.frame(k1, Matr2, SE, CILeft, CIRight)
cat(varlist[i], N, file = "Out.csv", append = TRUE)
cat("\n", file = "Out.csv", append = TRUE)
write.table(df, "Out.csv", sep = ",", col.names = !file.exists("Out.csv"), append = T)

#pdf(paste0('Figure',varlist[i],'.pdf'))
#pll <-   plot(CILeft, ylim = c(0,10),  type="l",
#     lwd =2,
#     col = "red", 
#     ylab = "Tail index estimate with confidence intervals",
#     xlab = "k",
#     main = paste0('Log-log rank-size plot ', varlist[i]), lty=2)
#ytick<-k1
#axis(side=2, at=ytick, labels = TRUE)
#points(CIRight, type="l", col = "red", lty=2)
#lines(Matr2, type="l", col = "blue")
#dev.off() 

tmp <- as.data.frame(df)
ggsave(filename=paste0('FigureInfHill',varlist[i],'.pdf'), width = 10, height = 7, scale = 1, 
       ggplot() +
         geom_line(tmp, mapping = aes(x = k1, y = Matr2), lwd =1,color="blue",) +
         geom_line(tmp, mapping = aes(x = k1, y = CILeft), lwd =1,color="red",linetype="dashed")+
         geom_line(tmp, mapping = aes(x = k1, y = CIRight), lwd =1,color="red",linetype="dashed")+
         ggtitle(paste0('Log-log rank-size plot ', varlist[i]))+
         ggeasy::easy_center_title()+coord_cartesian(xlim = c(2.5, 15), expand = 0)+
         labs(y= "Tail index estimate with confidence intervals", x = "k")+ 
         geom_point() + scale_x_continuous(n.breaks = 6))
#dev.copy(pll, filename=paste0('Figure',varlist[i-1],'.pdf'))
}

tmp <- as.data.frame(df)

ggplot() +
  geom_line(tmp, mapping = aes(x = k1, y = Matr2), lwd =1,color="blue",) +
  geom_line(tmp, mapping = aes(x = k1, y = CILeft), lwd =1,color="red",linetype="dashed")+
  geom_line(tmp, mapping = aes(x = k1, y = CIRight), lwd =1,color="red",linetype="dashed")+
  ggtitle(paste0('Log-log rank-size plot ', varlist[i]))+
  ggeasy::easy_center_title()+
  labs(y= "Tail index estimate with confidence intervals", x = "k")+ geom_point()+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))

ggplot() +
  geom_line(tmp, mapping = aes(x = k1, y = Matr2), lwd =1,color="blue",) +
  geom_line(tmp, mapping = aes(x = k1, y = CILeft), lwd =1,color="red",linetype="dashed")+
  geom_line(tmp, mapping = aes(x = k1, y = CIRight), lwd =1,color="red",linetype="dashed")+
  ggtitle(paste0('Log-log rank-size plot ', varlist[i]))+
  ggeasy::easy_center_title()+coord_cartesian(xlim = c(2.5, 15), expand = 0)+
  labs(y= "Tail index estimate with confidence intervals", x = "k")+ 
  geom_point() + scale_x_continuous(n.breaks = 6)


                     