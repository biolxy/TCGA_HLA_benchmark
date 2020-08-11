setwd('D:/project/HLA-benchmark/paper/manuscripts/FIGURE/revise/get_HR/20200412/univariable/p_05')
data <- read.table('notLOH_vs_LOH_absent.HR.plot.txt', header = TRUE, sep = "\t",stringsAsFactors=FALSE)
library(forestplot)

ci <- paste(data$CI_Low,"-",data$CI_High ,sep="")
pos_pts <- paste(data$pos,"(", data$pos_death, ")",sep="")
neg_pts <- paste(data$neg,"(", data$neg_death, ")",sep="")
## The rest of the columns in the table. 
tabletext <- cbind(c("Supertype",data$variable), 
                    #c("Cancer",data$Cancer),c("Stage",data$stage),
					#c("(-)allele.Pts.Num",data$other_allele.Pts),
					#c("(+)allele.Pts.Num",data$Allele.Pts),
					c("Condition", data$Condition),
					c("notLOH.Pts.Num(death)",pos_pts),
					c("Abs/LOH.Pts.Num(death)",neg_pts),
                    c("HR",data$HR), 
					c("95% CI",ci),				
                    c("P value",data$pvalue))
#png("Forestplot.png",width=1200, height=500)					

#png("Forestplot.png",width=1200, height=1200)
#pdf("Forestplot.pdf")
forestplot(labeltext=tabletext, graph.pos=5, 
           mean=c(NA,data$HR), 
           lower=c(NA,data$CI_Low), upper=c(NA,data$CI_High),
           xlab="<----Better overall survival----HR----Worse overall survival---->",
		   #xlab="Hazard Ratio",
		   txt_gp=fpTxtGp(label=gpar(cex=1.2),
                              ticks=gpar(cex=1.2),
                              xlab=gpar(cex =1.1),
                              title=gpar(cex = 1.2)),
           hrzl_lines= TRUE,
		   graphwidth = "auto", 
		   colgap = unit(1.8,"mm"),
		   lineheight = 'auto', 
		   col=fpColors(box="black", lines="black", zero = "gray50"),
		   zero=1, cex=0.12, boxsize=0.12,
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.12)