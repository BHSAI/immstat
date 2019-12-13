library("ggbiplot")

filetag = "BCF"

info = read.table("./data/info.csv", header = T, sep =';')
results = read.table(paste("./data/", filetag, "_output.txt", sep=''), header = T, sep =';')

d0 = merge(info, results, by="Subject_ID")	
d1 = subset(d0, Cohort != "Control")

print(length(d1$Subject_ID))

d_dfd = subset(d1, Cohort == "RRRx")
d_std = subset(d1, Cohort == "RRR")
d_new = rbind(d_dfd, d_std)

#figure options
plot_pca = function(d_new){

	ir.cohort = d_new$Cohort
#	imm_val = log(d_new[,5:12])
	
	Mem = d_new$Mem.PF.CSP.T6
	Naive = d_new$Naive.PF.CSP.T6
	ETL = d_new$ETL.PF.CSP.T6
	RM = d_new$RM.PF.CSP.T6
	AM = d_new$AM.PF.CSP.T6
	Ki67 = d_new$Ki67.PF.CSP.T6
	SM = d_new$SM.PF.CSP.T6
	UM = d_new$UM.PF.CSP.T6
	ET = d_new$ET.PF.CSP.T6
	LT = d_new$LT.PF.CSP.T6
	CD20 = d_new$CD20.PF.CSP.T6
	
	imm_val = data.frame(Mem, Naive, ETL, RM, AM, Ki67, SM, UM, ET, LT, CD20)
	
	#fill in missing values
	for (i in 1:length(imm_val)){	#columns
		for (j in 1:length(imm_val[,1])){	#subjects
			if (is.na(imm_val[j,i])){
				imm_val[j,i] = median(imm_val[,i], na.rm = TRUE)
				}
			}
		}
	

	ir.pca = prcomp(imm_val, center = TRUE, scale = TRUE)
	print(summary(ir.pca))
	dot_col = c("red", "blue")
	#biplot(ir.pca, scale = 0.5)
	#plot(ir.pca, type ='l')
	
	g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.cohort, ellipse = TRUE, 
              circle = TRUE)
	g <- g + scale_color_discrete(name = '')
	g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
	
}

cols = colnames(d_new)
plot_pca(d_new)
