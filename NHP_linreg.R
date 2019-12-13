source("immstat.R")

filetag = c("SERO", "FLR", "MESO", "TCF")
assayID = filetag
pathtag = "./results_adj/"

#extract columns to exclude
info = read.table("./data/info.csv", header = T, sep=',')
colnames_list = strsplit(colnames(info), " ")
exclude_col = array("null", dim = c(length(colnames_list)))
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

	
d_grp1 = subset(d_new, Adjuvant == "ALFA")
d_grp2 = subset(d_new, Adjuvant == "ALFQ")
d_grp3 = subset(d_new, Adjuvant == "ALFQA")	

n_col = length(exclude_col)
n_tot = length(colnames(d_new))
ALFA = array(0.0, dim = c(n_tot - n_col))
ALFQ = array(0.0, dim = c(n_tot - n_col))
ALFQA = array(0.0, dim = c(n_tot - n_col))
param = array("test", dim = c(n_tot - n_col))

for(i in 1:(n_tot - n_col)){
	param[i] = colnames(d_new)[i+n_col]
	ALFA[i] = median(d_grp1[,i+n_col], na.rm = TRUE)
	ALFQ[i] = median(d_grp2[,i+n_col], na.rm = TRUE)
	ALFQA[i] = median(d_grp3[,i+n_col], na.rm = TRUE)
	}
	
reg_data = data.frame(param, ALFA, ALFQ, ALFQA)
 
fit_adj = lm(ALFQA ~ ALFA + ALFQ + 0, data = reg_data)

summary(fit_adj)

write.table(reg_data, paste(pathtag, "linreg.csv", sep=''), sep=';', row.names = F)
	