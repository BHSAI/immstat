library("corrplot")
source("addicted.R")

# extract_timepoint = function(colname){
	# tp = 0
	# time = c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7")
	# for (i in 1:length(time)){
		# if (grepl(time[i], colname) == TRUE){
			# tp = time[i]
			# }
		# }
	# tp
	# }

extract_param = function(param_name){					#removes just timestamp
	spl = strsplit(param_name, ".", fixed = TRUE)[[1]]
	param = spl[1]
	for(i in 2:(length(spl)-1)){
		param = paste(param, spl[i], sep=".")
		}
	param
	}
	
extract_param2 = function(param_name){					#removes timestamp and assayID
	spl = strsplit(param_name, ".", fixed = TRUE)[[1]]
	param = spl[2]
	for(i in 3:(length(spl)-1)){
		param = paste(param, spl[i], sep=".")
		}
	param
	}
	
add_assayID = function(d1, assay_ID){
	exclude_col = c("Subject_ID", "Cohort", "Status")
	new_cols = colnames(d1)
	for(i in 1:length(new_cols)){
		skip = FALSE
		if (match(new_cols[i], exclude_col, nomatch = 0) != 0){
						skip = TRUE
						}
		if (skip == FALSE){
			new_cols[i]=paste(assay_ID,new_cols[i], sep='.')
			}
		}
	colnames(d1) = new_cols
	d1
	}
	
vacc_induced = function(param, q1){
	pass = FALSE
	
	#grab all values with param
	i = match(extract_param(param), q1$cols)

	time_point = c("T2", "T4", "T5", "T6", "T7")
	
	class_sum = 0
	for(j in time_point){
		p_var = paste(j,"vacc.class", sep='.')
		class_sum = class_sum + q1[[p_var]][i]
		print
		print(q1[[p_var]][i])
		}
	
	if (class_sum > 0){
		pass = TRUE
		}
	pass
	}
	
group_diff = function(param, q1){
	pass = FALSE
	i = match(param, q1$cols)
	if (is.na(i) == FALSE){
		p = q1$pval[i]
		q = q1$q_val[i]

		if(is.na(p) == FALSE && is.na(q) == FALSE){
			if(vacc_induced(param, q1) == TRUE){
				if (min(p)<0.05 && min(q)< 0.20){
					pass = TRUE
					}
				}
			}
		}
	pass
	}
	
	
plot_matrix = function(mcor, filename){
	varnum = dim(mcor)[1]
	fig_width = 5+round(varnum/8)
	png(filename, width = fig_width, height = fig_width-1, units = 'in', res = 300)
	corrplot(mcor, type="upper", order="hclust", tl.col="black", tl.srt=45)
	#heatmap(mcor, scale='none')
	dev.off()
	}
	
generate_clusters = function(mcor){
	dissimilarity = 1- abs(mcor)
	distance = as.dist(dissimilarity)
	hc = hclust(distance)
	hc
	}
	
plot_dendrogram = function(hc, filename, clust_ID){
	#generate dendrogram
	#dend = as.dendrogram(hc)
	varnum = length(clust_ID)
	
	fig_width = 5+round(varnum/8)
	png(filename, width = fig_width, height = 7, units = 'in', res = 300)
	#plot(hc)
	cluster_num = max(clust_ID)
	color_array = c("black", "red", "green", "blue", "orange", "purple")
	for(i in 1:10){
		color_array = c(color_array, color_array)
		}
	A2Rplot(hc, k = cluster_num, boxes = FALSE, col.up = "gray50", col.down = color_array[1:cluster_num])
	dev.off()
}

corr_analysis = function(results_table){

	exclude_col = c("Subject_ID", "Cohort", "Status", "vacc")
	d1=subset(results_table, Cohort != "Control")

	Subject_ID = d1$Subject_ID
	
	new_data = data.frame(Subject_ID)
	cols = colnames(d1)
	col_num = dim(d1)[2]

	for(i in 1:col_num){
		print(cols[i])
		skip = F
		if (match(cols[i], exclude_col, nomatch = 0) != 0){ #skip certain columns like SubjectID
			skip = T
			}
		if (skip == F){
			new_name = cols[i]
			new_data[[new_name]] = d1[,i]
			}
		}
		
	output_table = "NA"
			
	new_data$Subject_ID = NULL
	if (length(new_data) > 0){
		mcor = cor(new_data, use = "pairwise.complete.obs", method = "spearman")
		col_num = dim(new_data)[2]
		threshold = 0.05
		for(i in 1:col_num){
			for(j in 1:col_num){
				p = cor.test(new_data[,i], new_data[,j])$p.value
				if (p > 0.05){
					mcor[i,j] = 0
					}
				}
			}
		}
		
	mcor
	}

#build data table for analysis
info = read.table("./data/info.csv", header = T, sep =';')
BCF = read.table("./data/BCF_output.txt", header = T, sep =';')
IGG = read.table("./data/IGG_output.txt", header = T, sep =';')
ICC = read.table("./data/ICC_output.txt", header = T, sep =';')
BELI = read.table("./data/BELI_output.txt", header = T, sep =';')

#AELI = read.table("./data/AELI_output.txt", header = T, sep =';')
#CYT = read.table("./data/CYT_output.txt", header = T, sep =';')
#SERO = read.table("./data/SERO_output.txt", header = T, sep =';')

exclude_col = c("Subject_ID", "Cohort", "Status")

#since different assays are being combined, add assay tag
BCF = add_assayID(BCF, "BCF")
IGG = add_assayID(IGG, "IGG")
ICC = add_assayID(ICC, "ICC")
BELI = add_assayID(BELI, "BELI")

#AELI = add_assayID(AELI, "AELI")
#CYT = add_assayID(CYT, "CYT")
#SERO = add_assayID(SERO, "SERO")

#combine data sets
d0 = merge(info, BCF, by="Subject_ID")	
d1 = merge(d0, IGG, by="Subject_ID")
d_final = merge(d1, BELI, by="Subject_ID")	
#d_final = merge(d1, ICC, by="Subject_ID")
 

#d3 = merge(d2, AELI, by="Subject_ID")
#d_final = merge(d4, CYT, by="Subject_ID")
#d_final = merge(d5, SERO, by="Subject_ID")

#d_final = merge(info, ICC, by="Subject_ID")

#compile p and q values from merged data
#filetag = c("ICC", "IGG", "BELI", "BCF")
#filetag = c("ICC")
filetag = c("IGG", "BELI", "BCF")
pathtag = "./results_group_qval/"

q0 = read.table(paste(pathtag,filetag[1],".merged.csv", sep=''), sep=',', header = T)
q1 = q0
if(length(filetag)>1){
	for(i in 2:length(filetag)){
		q0 = read.table(paste(pathtag,filetag[i],".merged.csv", sep=''), sep=',', header = T)
		q1 = rbind(q0, q1)
		}
	}

#filter data by vaccine induced and timepoint
cols = colnames(d_final)
col_num = dim(d_final)[2]

#analysis parameters
time_val = "T5"
split_time = TRUE
plot_dendro = TRUE
plot_matrix = FALSE

for(i in 4:col_num){
	if(vacc_induced(cols[i], q1) == FALSE){
		d_final[[cols[i]]] <- NULL
		}
	if(split_time == TRUE){
			if (grepl(time_val, cols[i]) == FALSE){
				d_final[[cols[i]]] <- NULL
				}
			}
}

#clean param names
cols = colnames(d_final)
for(i in cols){
	skip = FALSE
	exclude_col = c("Subject_ID", "Cohort", "Status")
	if (match(i, exclude_col, nomatch = 0) != 0){
					skip = TRUE
					}
	if(skip == FALSE){
		new_name = extract_param(i)
		d_final[[new_name]] = d_final[[i]]
		d_final[[i]] = NULL
		}
	}
	
# d_final$cols = NULL
# d_final$cols = d_final$new_name
# d_final$new_name = NULL

#calculate correlation matrix
mcor = corr_analysis(d_final)

#output results
if(plot_matrix == TRUE){
	plot_matrix(mcor, paste("matrix",time_val,"png", sep ="."))
	}

#delineate clusters
hc = generate_clusters(mcor)
clust_ID = cutree(hc, h = 0.5)
cols = labels(clust_ID)

#plot dendrogram
if(plot_dendro == TRUE){
	plot_dendrogram(hc, paste("dendro",time_val,"png", sep ="."), clust_ID)
	}

#output cluster IDs
write.table(data.frame(cols,clust_ID), file = paste("clustID",time_val,"txt", sep ="."), append=FALSE, sep =';', row.names=FALSE)

#summary p-values for manuscript
fig_table = read.table("dendro_xname.csv", header = TRUE, sep=';')

#make table with param name, vacc p+q, group p+q, prot p+q for T0, T2, T4, T5, T6, T7
time_point = c("T2", "T4", "T5", "T6", "T7")
type = c("vacc", "group", "prot")
n_param = length(fig_table$x)

for(i in 1:n_param){
	#param_name = extract_param(as.character(fig_table$x[i]))
	param_name = as.character(fig_table$x[i])
	fig_table$param[i] = param_name
	q2 = subset(q1, cols == param_name) #grab param from data 
	for(j in type){
		for(k in time_point){				#extract timepoint specific data 
			p_class = paste(k, j, "class", sep='.')
			p_var = paste(k, j, "p", sep='.')
			
			fig_table[[p_var]][i] = 1.0
			if(q2[[p_class]][1] == 1){
				fig_table[[p_var]][i] = q2[[p_var]][1]
				}	
			}
		}
	}
	fig_table$x = NULL
	
#draw figure
filename = "fig_vacc.png"
png(filename, width = 10, height = 5, units = 'in', res = 300)

plot(1, type="n", xlab ="", ylab = "Timepoint", xlim = c(1, n_param), ylim = c(0,22), axes=FALSE)
for(i in 1:n_param){
	for(j in 1:length(time_point)){
		xval = i
		yval = 21-j
		text(-1, yval, labels=time_point[j], cex=0.5)
		text(n_param+2, yval, labels=time_point[j], cex=0.5)
		abline(h=yval-0.5, col = 1)
		abline(h=yval+0.5, col = 1)

		#vaccine induced responses
		pt_col = "white"
		p_var = paste(time_point[j], "vacc.p", sep='.')

		if(fig_table[[p_var]][i] < 0.05){
			pt_col = "gray85"
			} 
		if(fig_table[[p_var]][i] < 0.01){
			pt_col = "gray50"
			} 
		if(fig_table[[p_var]][i] < 0.001){
			pt_col = "black"
			} 
			
		if(pt_col != "white"){	
			points(xval, yval, type='o', col= pt_col, pch = 16, cex = 0.75)
			}
		
		#group level differences
		yval = 14-j
		text(-1, yval, labels=time_point[j], cex=0.5)
		text(n_param+2, yval, labels=time_point[j], cex=0.5)
		abline(h=yval-0.5, col = 1)
		abline(h=yval+0.5, col = 1)
		
		pt_col = "white"
		p_var = paste(time_point[j], "group.p", sep='.')

		if(fig_table[[p_var]][i] < 0.05){
			pt_col = "cornflowerblue"
			} 
		if(fig_table[[p_var]][i] < 0.01){
			pt_col = "blue"
			} 
		if(fig_table[[p_var]][i] < 0.001){
			pt_col = "blue4"
			} 
			
		if(pt_col != "white"){	
			points(xval, yval, type='o', col= pt_col, pch = 16, cex = 0.75)
			}
			
		#protection level differences
		yval = 7-j
		text(-1, yval, labels=time_point[j], cex=0.5)
		text(n_param+2, yval, labels=time_point[j], cex=0.5)
		abline(h=yval-0.5, col = 1)
		abline(h=yval+0.5, col = 1)
		
		pt_col = "white"
		p_var = paste(time_point[j], "prot.p", sep='.')

		if(fig_table[[p_var]][i] < 0.05){
			pt_col = "indianred"
			} 
		if(fig_table[[p_var]][i] < 0.01){
			pt_col = "red"
			} 
		if(fig_table[[p_var]][i] < 0.001){
			pt_col = "red4"
			} 
			
		if(pt_col != "white"){	
			points(xval, yval, type='o', col= pt_col, pch = 16, cex = 0.75)
			}
		}
	}
	
dev.off()



