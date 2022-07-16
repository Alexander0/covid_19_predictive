names(LogLogWork)
varlist <- names(LogLogWork)[2:25]

for (i in 2:(length(varlist)+1)){
  
N<-sum(LogLogWork[!is.na(LogLogWork[,i]),i]>0)
N
Matr2<-replicate(N-4, 0)
SE<-replicate(N-4, 0)
CILeft<-replicate(N-4, 0)
CIRight<-replicate(N-4, 0)
for (k in 5:N){
  lin1<-lm(LogLogWork$lRank[1:k]~as.matrix(LogLogWork[1:k,i]), data=LogLogWork)
  C<-coef(lin1)
  Matr2[k-4]<- -C[[2]]
  SE[k-4]<- -C[[2]]*sqrt(2)/sqrt(k)
  CILeft[k-4]<- -C[[2]]+1.96*C[[2]]*sqrt(2)/sqrt(k)
  CIRight[k-4]<- -C[[2]]-1.96*C[[2]]*sqrt(2)/sqrt(k)
}

k1<-5:N
Matr2  
SE
CILeft
CIRight
df <- data.frame(k1, Matr2, SE, CILeft, CIRight)
cat("Austray", N, file = "Out.csv", append = TRUE)
cat("\n", file = "Out.csv", append = TRUE)
write.table(df, "Out.csv", sep = ",", col.names = !file.exists("Out.csv"), append = T)

plot(CILeft, ylim = c(0,4),  type="l",
     lwd =2,
     col = "red", 
     ylab = "Tail index estimate with confidence intervals",
     xlab = "k",
     main = "Log-log rank-size plot - Australiay", lty=2)
ytick<-seq(0, 30, by=1)
axis(side=2, at=ytick, labels = FALSE)
points(CIRight, type="l", col = "red", lty=2)
lines(Matr2, type="l", col = "blue")
}
