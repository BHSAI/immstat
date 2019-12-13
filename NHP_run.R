library(beeswarm)
library("corrplot")
source("addicted.R")
source("immstat.R")
library("ggbiplot")
library("caret")


filetag = c("SERO", "FLR", "MESO", "TCF")
assayID = filetag
skipT0 = c(FALSE, FALSE, TRUE, TRUE)
draw_fig = FALSE
pathtag = "./results_adj/"

#extract columns to exclude
info = read.table("./data/info.csv", header = T, sep=',')
colnames_list = strsplit(colnames(info), " ")
exclude_col = array("null", dim = c(length(colnames_list)))
for(i in 1:length(exclude_col)){
	exclude_col[i] = colnames_list[i][[1]]
	}
	
#DESCRIPTIVE STATISTICS
for (i in 1:length(filetag)){
	print(paste("RUNNING ANALYSIS FOR ", filetag[i], sep=''))
	pvals = run_analysis(filetag[i], pathtag, draw_fig, skipT0[i], assayID[i], exclude_col)
	calc_qvalue(pvals, filetag[i], pathtag, exclude_col)
	}

#CORRELATION ANALYSIS
print("CORRELATION ANALYSIS")
d0 = read.table(paste("./data/info.csv", sep=''), header = T, sep =',')
for(i in filetag){
	d1=read.table(paste("./data/", i, ".csv", sep=''), header = T, sep =',')
	d1 = add_assayID(d1, i, exclude_col)
	d0 = merge(d0, d1, by="Subject_ID")
	}

#calculate correlation matrix
d_final = d0
exclude_time = "T0"
mcor = corr_analysis(d_final, exclude_col, exclude_time)					
plot_matrix(mcor,  paste(pathtag, "matrix.png", sep =''))

#delineate clusters
hc = generate_clusters(mcor)
clust_ID = cutree(hc, h = 0.6)
cols = labels(clust_ID)

#plot dendrogram
plot_dendrogram(hc, paste(pathtag, "dendro.png", sep =''), clust_ID)

#output cluster IDs and merge with statistics
clust_list = data.frame(cols, clust_ID)

#compile statistics
c0 = read.table(paste(pathtag, filetag[1], ".qvalue.csv", sep=''), header = T, sep =';')

for(i in 2:length(filetag)){
	c1 = read.table(paste(pathtag, filetag[i], ".qvalue.csv", sep=''), header = T, sep =';')
	c0 = rbind(c1, c0)
	}

c1 = subset(c0, timepoint != "T0")
c_output = merge(c1, clust_list, by="cols")
write.table(c_output, paste(pathtag, "all_data.csv", sep=''), sep=';', row.names = F)

#draw dendro statistics
draw_dendro_stats_vacc(pathtag, c_output)
draw_dendro_stats_group(pathtag, c_output)

#PRINCIPAL COMPONENT ANALYSIS
print("PRINCIPAL COMPONENT ANALYSIS")
Adjuvant = d_final$Adjuvant
Subject_ID = d_final$Subject_ID
TCF.LN.CD4 = d_final$TCF.ln.AgSp.CD4.CSP.T5
Meso.IL4 = d_final$MESO.liver.IL4.CSP.T5
Meso.IL5 = d_final$MESO.liver.IL5.CSP.T5
Meso.Multi = d_final$MESO.liver.IFNy.CSP.T5
Meso.PBMC.IL4 = d_final$MESO.pbmc.IL4.FMP14.T5
FLR.IL2 = d_final$FLR.IL2.CSP.T5
FLR.IFNy = d_final$FLR.IFNy.CSP.T5
FLR.IL2.T3 = d_final$FLR.IL2.CSP.T3
FLR.IFNy.T3 = d_final$FLR.IFNy.CSP.T3
SERO.ELISA = d_final$SERO.ELISA.NANP.T5

imm_val = data.frame(TCF.LN.CD4, Meso.IL4, Meso.IL5, Meso.Multi, Meso.PBMC.IL4, FLR.IL2, FLR.IFNy,FLR.IL2.T3, FLR.IFNy.T3, SERO.ELISA)

group = array("Saline+None", dim = c(length(d_final$Cohort)))
for(i in 1:length(d_final$Cohort)){
	group[i] = paste(d_final$Cohort[i], d_final$Adjuvant[i], sep='+')
	}
	
dot_col = c("red", "blue", "green", "yellow", "purple", "pink", "orange", "black") #based on number of cohorts to color
plot_pca(pathtag, "test_pca_adj.png", imm_val, d_final$Adjuvant, dot_col)

# RANDOM FOREST MODEL
print("RANDOM FOREST MODEL")
#c_output is param-indexed data
#d_final is subject-indexed data

#subset of data 
#d_rf = data.frame(Subject_ID, Adjuvant, TCF.LN.CD4, Meso.IL4, Meso.IL5, Meso.Multi, Meso.PBMC.IL4, FLR.IL2, FLR.IFNy,FLR.IL2.T3, FLR.IFNy.T3, SERO.ELISA)
d_grp1 = subset(d_final, Adjuvant == "ALFA")
d_grp2 = subset(d_final, Adjuvant == "ALFQ")

# d_grp1 = subset(d_rf, Adjuvant == "ALFA")
# d_grp2 = subset(d_rf, Adjuvant == "ALFQ")
d_new = rbind(d_grp1, d_grp2)
d_new$Adjuvant = factor(d_new$Adjuvant)

#keep only vaccine-induced responses
cols = colnames(d_new)
col_num = dim(d_new)[2]	
for(i in 1:col_num){
	if(exclude_col_filter(cols[i], exclude_col) == FALSE){ 
		if(group_diff(cols[i], c_output) == TRUE && vacc_induced(cols[i], c_output) == TRUE){
		#if(vacc_induced(cols[i], c_output) == TRUE){
			print(cols[i])
			} else {
				d_new[[cols[i]]] <- NULL
			}
		}
	}
	
for(j in exclude_col){
		if(j != "Adjuvant" && j != "Subject_ID"){
			d_new[[j]] <- NULL
			}
		}	

#fill in missing values
var_names = colnames(d_new)
for (i in 1:length(var_names)){	#columns
	if(exclude_col_filter(var_names[i], exclude_col) == FALSE){
		for (j in 1:length(d_new[,1])){	#subjects
			if (is.na(d_new[j,i])){
				d_new[j,i] = median(d_new[,i], na.rm = TRUE)
				}
			}
		}
	}

#center and scale
var_names = colnames(d_new)
for(i in 1:length(var_names)){
	if(exclude_col_filter(var_names[i], exclude_col) == FALSE){
		if(sum(d_new[[var_names[i]]])>0){
			d_new[[var_names[i]]] = scale(d_new[[var_names[i]]], center=TRUE, scale=TRUE)
			}
		}
	}

#run Random Forest Model
ctrl = trainControl(method="repeatedcv", number=5, repeats=10, selectionFunction = "oneSE")

#leave one out analysis
id = unique(d_new$Subject_ID)
correct = 0
prediction_table = d_new
prediction_table$pred.adjuvant = d_new$Adjuvant #give it the same levels/format

print("leave-one-out analysis")
for (i in 1:length(id)){
	train_data = d_new[d_new$Subject_ID!=id[i],]
	test_data = d_new[d_new$Subject_ID==id[i],]
	train_data$Subject_ID <- NULL
	test_data$Subject_ID <- NULL

	#generate random forest with training data
	# trf = train(Adjuvant ~ .,
	# data=train_data, method="rf", metric="Kappa", 
	# trControl=ctrl)
	
	trf = train(Adjuvant ~ .,	
		data=train_data, method="ranger", metric="Kappa", 
	trControl=ctrl, importance = "permutation")
	
	#predict on leave-one-out subject
	prediction_table$pred.adjuvant[i] = predict(trf, test_data, "raw")
	}

confusion = confusionMatrix(prediction_table$pred.adjuvant, prediction_table$Adjuvant)
print(confusion)

#output weights
print("generating tree on whole data set")
test_data_all = d_new
test_data_all$Subject_ID <- NULL
			
# trf = train(Adjuvant ~ .,
	# data=test_data_all, method="rf", metric="Kappa",
		# trControl=ctrl)
		
trf = train(Adjuvant ~ .,	
		data=train_data, method="ranger", metric="Kappa", 
	trControl=ctrl, importance = "permutation")

#variable importance
print(varImp(trf))
