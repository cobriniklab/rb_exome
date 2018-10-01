#!/usr/bin/Rscript


#filepaths


library(ggplot2)

somePDFPath = "C:\\temp\\some.pdf"
pdf(file=somePDFPath)  

for (i in seq(5,10))   
{   
  par(mfrow = c(2,1))
  VAR1=rnorm(i)  
  VAR2=rnorm(i)  
  plot(VAR1,VAR2)   
} 
dev.off() 
