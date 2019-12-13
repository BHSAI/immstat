source("immstat.R")

#extract columns to exclude
info = read.table("./data/info.csv", header = T, sep=',')
colnames_list = strsplit(colnames(info), " ")
exclude_col = array("null", dim = c(length(colnames_list)))

#define resample function
resample_data = function(a_grp1){
	d_grp1 = a_grp1
	n_col = length(exclude_col)
	n_tot = length(colnames(a_grp1))
	
	nsub = length(a_grp1$Adjuvant)
	for(i in 1:(n_tot- n_col)){
		for(j in 1:nsub){
			d_grp1[j, i+n_col] = a_grp1[sample(1:nsub,1),i+n_col]
			}
		}
	d_grp1
	}

filetag = c("SERO", "FLR", "MESO", "TCF")
assayID = filetag
pathtag = "./results_adj/"

for(i in 1:length(exclude_col)){
	exclude_col[i] = colnames_list[i][[1]]
	}

d0 = read.table(paste("./data/info.csv", sep=''), header = T, sep =',')
for(i in filetag){
	d1=read.table(paste("./data/", i, ".csv", sep=''), header = T, sep =',')
	d1 = add_assayID(d1, i, exclude_col)
	d0 = merge(d0, d1, by="Subject_ID")
	}

#data indexed by subject
d_final = d0

#param index data
c_output = read.table(paste(pathtag, "all_data.csv", sep=''), sep=';', header = T)

#keep only vaccine-induced and adj-different responses
d_new = d_final
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
			d_new[[var_names[i]]] = scale(d_new[[var_names[i]]], center=FALSE, scale=TRUE)
			}
		}
	}
	
a_grp1 = subset(d_new, Adjuvant == "ALFA")
a_grp2 = subset(d_new, Adjuvant == "ALFQ")
a_grp3 = subset(d_new, Adjuvant == "ALFQA")	

#resample data
sample_rounds = 10

for(i in 1:(n_tot - n_col)){
		param[i] = colnames(d_new)[i+n_col]
		}


ALFA_all = data.frame(param)
ALFQ_all = data.frame(param)
ALFQA_all = data.frame(param)

coeff_ALFA = array(0.0, dim = c(sample_rounds))
coeff_ALFQ = array(0.0, dim = c(sample_rounds))

for(k in 1:sample_rounds){
	d_grp1 = resample_data(a_grp1)
	d_grp2 = resample_data(a_grp2)
	d_grp3 = resample_data(a_grp3)

	#generate median values
	n_col = length(exclude_col)
	n_tot = length(colnames(d_new))
	ALFA = array(0.0, dim = c(n_tot - n_col))
	ALFQ = array(0.0, dim = c(n_tot - n_col))
	ALFQA = array(0.0, dim = c(n_tot - n_col))
	param = array("test", dim = c(n_tot - n_col))

	#select median values as representative values for each parameter
	for(i in 1:(n_tot - n_col)){
		ALFA[i] = median(d_grp1[,i+n_col], na.rm = TRUE)
		ALFQ[i] = median(d_grp2[,i+n_col], na.rm = TRUE)
		ALFQA[i] = median(d_grp3[,i+n_col], na.rm = TRUE)
		}
		
	reg.data = data.frame(ALFA, ALFQ, ALFQA)
	fit_adj = lm(ALFQA ~ ALFA + ALFQ + 0, data = reg_data)
	print(fit_adj$coefficients)
	coeff_ALFA[k] = as.numeric(fit_adj$coefficients[1])
	coeff_ALFQ[k] = as.numeric(fit_adj$coefficients[2])
		
	ALFA_all[[paste("ALFA",k,sep='.')]] = ALFA
	ALFQ_all[[paste("ALFQ",k,sep='.')]] = ALFQ
	ALFQA_all[[paste("ALFQA",k,sep='.')]] = ALFQA
	}
	
for(i in 1:length(param)){

	ALFA_array = as.numeric(ALFA_all[i,2:(sample_rounds+1)])
	ALFQ_array = as.numeric(ALFQ_all[i,2:(sample_rounds+1)])
	ALFQA_array = as.numeric(ALFQA_all[i,2:(sample_rounds+1)])

	ALFA_all$average[i] = mean(ALFA_array)
	ALFQ_all$average[i] = mean(ALFQ_array)
	ALFQA_all$average[i] = mean(ALFQA_array)
	ALFA_all$perc25[i] = quantile(ALFA_array)[2]
	ALFQ_all$perc25[i] = quantile(ALFQ_array)[2]
	ALFQA_all$perc25[i] = quantile(ALFQA_array)[2]
	ALFA_all$perc75[i] = quantile(ALFA_array)[4]
	ALFQ_all$perc75[i] = quantile(ALFQ_array)[4]
	ALFQA_all$perc75[i] = quantile(ALFQA_array)[4]
	}
	
linreg_data = data.frame(coeff_ALFA, coeff_ALFQ)

write.table(linreg_data, paste(pathtag, "linreg_sample.csv", sep=''), sep=';', row.names = F)
write.table(ALFA_all, paste(pathtag, "ALFA_sample.csv", sep=''), sep=';', row.names = F)
write.table(ALFQ_all, paste(pathtag, "ALFQ_sample.csv", sep=''), sep=';', row.names = F)
write.table(ALFQA_all, paste(pathtag, "ALFQA_sample.csv", sep=''), sep=';', row.names = F)

#write.table(reg_data, paste(pathtag, "linreg.csv", sep=''), sep=';', row.names = F)