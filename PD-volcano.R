data<-read.csv('1.csv',stringsAsFactors = F)
gene.list<-data$Symbol
fc.list<-data$es
p.value.list<-data$p
plot(fc.list,-log10(p.value.list),pch=8,xlab = 'effect size',ylab = '-log10 (p.value)',type = 'n')

for(i in c(1:length(gene.list))){
  if(p.value.list[i]<0.05){
    if(fc.list[i]>0.2){
      point.col<-'red'
    }else if(fc.list[i]<(-0.2)){
      point.col<-'blue'
    }else{
      point.col<-'grey'
    }
  }else{
    point.col<-'grey'
  }
  points(fc.list[i],-log10(p.value.list[i]),pch=20,col=point.col,cex=0.5)
}