library(beeswarm)
library("ggbiplot")

#figure options
# class0 = "Saline"
# class1 = "10ug FMP014"	
# class2 = "40ug FMP014"

class0 = "None"
class1 = "ALFA"
class2 = "ALFQ"

time_array = c("T0", "T1", "T2", "T3", "T4", "T5")

extract_timepoint = function(colname){
	tp = 0
	for (i in 1:length(time_array)){
		if (grepl(time_array[i], colname) == TRUE){
			tp = time_array[i]
			}
		}
	tp
	}

extract_assayID = function(param_name){
	spl = strsplit(param_name, ".", fixed = TRUE)[[1]]
	param = spl[1]
	param
	}
	
calc_effect_size = function(grp1, grp2, dat){
	mean1 = mean(grp1, na.rm = TRUE)
	mean2 = mean(grp2, na.rm = TRUE)
	mean_sd = mean(c(sd(grp1, na.rm = TRUE), sd(grp2, na.rm = TRUE)))
	
	effect_size = abs(mean1 - mean2)/mean_sd
	effect_size
	}
	
add_assayID = function(d1, assay_ID, exclude_col){
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
	
vacc_induced = function(param, q1){
	pass = FALSE
	
	i = match(param, q1$cols)
	#param_data = subset(q1, cols == fig_table$param[i])
	if (is.na(i) == FALSE){
		p_vals = c(1.0, q1$vacc.all.p[i], q1$vacc.grp2.p[i],q1$vacc.grp3.p[i])
		q_vals = c(1.0, q1$vacc.all.q[i], q1$vacc.grp2.q[i],q1$vacc.grp3.q[i])

		if (min(p_vals, na.rm = TRUE)<0.05 && min(q_vals, na.rm = TRUE)< 0.20){
			pass = TRUE
			}
			
		}
	pass
	}
	
group_diff = function(param, q1){
	pass = FALSE
	
	i = match(param, q1$cols)
	if (is.na(i) == FALSE){
		p_vals = c(1.0, q1$group.p[i])
		q_vals = c(1.0, q1$group.q[i])

		if (min(p_vals, na.rm = TRUE)<0.05 && min(q_vals, na.rm = TRUE)< 0.20){
			pass = TRUE
			}
		}
	pass
	}
	
exclude_col_filter = function(param, exclude_col){	
	pass = FALSE
	if (match(param, exclude_col, nomatch = 0) != 0){
		pass = TRUE
		}
	pass
	}
	
# group_diff = function(param, q1){
	# pass = FALSE
	# i = match(param, q1$cols)
	# if (is.na(i) == FALSE){
		# p = q1$group.p[i]
		# q = q1$group.q[i]

		# if(is.na(p) == FALSE && is.na(q) == FALSE){
			# if(vacc_induced(param, q1) == TRUE){
				# if (min(p)<0.05 && min(q)< 0.20){
					# pass = TRUE
					# }
				# }
			# }
		# }
	# pass
	# }

#CORRELATION ANALYSIS FUNCTIONS
plot_matrix = function(mcor, filename){
	varnum = dim(mcor)[1]
	fig_width = 5+round(varnum/7)
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

draw_dendro_stats_vacc= function(filepath, merge_data){
	filename = paste(filepath, "dendro_stat_vacc.png", sep='')
	fig_table = read.table("dendro_xname.csv", header = TRUE, sep=';')
	n_param = length(fig_table$x)
	
	png(filename, width = 10, height = 2, units = 'in', res = 300)
	plot(1, type="n", xlab ="", ylab = "", xlim = c(1, n_param), ylim = c(0,3), axes=FALSE)
	
	fig_table$p = 1.0
	fig_table$q = 1.0
	
	for(i in 1:n_param){
		
		param_name = as.character(fig_table$x[i])
		fig_table$param[i] = param_name
		
		param_data = subset(merge_data, cols == fig_table$param[i])
		p_vals = c(1.0, param_data$vacc.all.p[1], param_data$vacc.grp2.p[1],param_data$vacc.grp3.p[1])#not being read as numbers
		q_vals = c(1.0, param_data$vacc.all.q[1], param_data$vacc.grp2.q[1],param_data$vacc.grp3.q[1])
		fig_table$p[i] = min(p_vals, na.rm= TRUE)
		fig_table$q[i] = min(q_vals, na.rm=TRUE)
		
		#draw on vacc induced stats
		pt_col = "white"
		if(fig_table$p[i] < 0.05 && fig_table$q[i] < 0.20){
			pt_col = "gray85"
			}
		if(fig_table$p[i] < 0.01 && fig_table$q[i] < 0.20){
			pt_col = "gray50"
			}
		if(fig_table$p[i] < 0.001 && fig_table$q[i] < 0.20){
			pt_col = "black"
			}
		fig_table$pt_col[i] = pt_col
		if(pt_col != "white"){	
			points(i, 2, type='o', col= pt_col, pch = 16, cex = 0.75)
			}
		write.table(fig_table, "debug.csv", sep =';', row.names = F)
			
		}
		
	dev.off()
	}
	
draw_dendro_stats_group= function(filepath, merge_data){
	filename = paste(filepath, "dendro_stat_group.png", sep='')
	fig_table = read.table("dendro_xname.csv", header = TRUE, sep=';')
	n_param = length(fig_table$x)
	
	png(filename, width = 10, height = 2, units = 'in', res = 300)
	plot(1, type="n", xlab ="", ylab = "", xlim = c(1, n_param), ylim = c(0,3), axes=FALSE)
	
	fig_table$p = 1.0
	fig_table$q = 1.0
	
	for(i in 1:n_param){
		
		param_name = as.character(fig_table$x[i])
		fig_table$param[i] = param_name
		
		param_data = subset(merge_data, cols == fig_table$param[i])
		p_vals = c(1.0, param_data$group.p[1])#not being read as numbers
		q_vals = c(1.0, param_data$group.q[1])
		fig_table$p[i] = min(p_vals, na.rm= TRUE)
		fig_table$q[i] = min(q_vals, na.rm=TRUE)
		
		#draw on vacc induced stats
		pt_col = "white"
		if(fig_table$p[i] < 0.05 && fig_table$q[i] < 0.20){
			pt_col = "cornflowerblue"
			}
		if(fig_table$p[i] < 0.01 && fig_table$q[i] < 0.20){
			pt_col = "blue"
			}
		if(fig_table$p[i] < 0.001 && fig_table$q[i] < 0.20){
			pt_col = "blue4"
			}
		fig_table$pt_col[i] = pt_col
		if(pt_col != "white"){	
			points(i, 2, type='o', col= pt_col, pch = 16, cex = 0.75)
			}
		write.table(fig_table, "debug.csv", sep =';', row.names = F)
			
		}
		
	dev.off()
	}


corr_analysis = function(results_table, exclude_col, exclude_time){

	d1=subset(results_table, Adjuvant != class0)

	Subject_ID = d1$Subject_ID
	
	new_data = data.frame(Subject_ID)
	cols = colnames(d1)
	col_num = dim(d1)[2]

	for(i in 1:col_num){
		skip = F
		if (match(cols[i], exclude_col, nomatch = 0) != 0){ #skip certain columns like SubjectID
			skip = T
			}
		if (match(extract_timepoint(cols[i]), exclude_time, nomatch =0) != 0){ #skip certain time points like T0
			skip = T
			}
		if (skip == F){
			new_name = cols[i]
			new_data[[new_name]] = d1[,i]
			}
		}
		
	output_table = "NA"
			
	new_data$Subject_ID = NULL
	debug_t=colnames(new_data)
	if (length(new_data) > 0){
		mcor = cor(new_data, use = "pairwise.complete.obs", method = "spearman")
		col_num = dim(new_data)[2]
		threshold = 0.05
		for(i in 1:col_num){
			for(j in 1:col_num){
				
				p = cor.test(new_data[,i], new_data[,j])$p.value
				if(is.na(p) == TRUE){
					mcor[i,j] = 0
					} 
				if (is.na(p) == FALSE && p > 0.05){
					mcor[i,j] = 0
					}
				}
			}
		}
		
	mcor
	}
#CORRELATION ANALYSIS FUNCTION END

#DESCRIPTIVE STATISTICS FUNCTIONS
plot_chart = function(dat, d1, filetag, pathtag, draw_fig, skipT0){

	#define groups
	# group1=subset(d1, Cohort == class0)
	# group2=subset(d1, Cohort == class1)
	# group3=subset(d1, Cohort == class2)
	
	group1=subset(d1, Adjuvant == class0)
	group2=subset(d1, Adjuvant == class1)
	group3=subset(d1, Adjuvant == class2)
	
	# a0 = rbind(group1, group2)
	# a1 = rbind(group3, a0)
	# all = a1
	all = d1
	dot_col = c("red", "blue", "green", "gray")

	plot_data = data.frame(all$Cohort)
	plot_data$dat = all[[dat]]
	plot_data$Cohort = factor(all$Cohort)
	plot_data$Adjuvant = factor(all$Adjuvant) 
	plot_data$Dose = factor(all$Dose)
	p = 1.0
	
	n1 = 0
	n2 = 0
	output_val = c(1.0, 0.0, 0.0,1.0,1.0,1.0,0.0, 0.0, 0.0)
	#output_val 1 p-value
	#output_val 2 magnitude
	#output_val 3 effect size
	#output_val 4 p-value group2 vs. T0
	#output_val 5 p-value group3 vs. T0
	#output_val 6 p-value all vs. T0
	#output_val 7 effect size group2 vs group1
	#output_val 8 effect size group3 vs group1
	#output_val 9 effect size all vs group1
	
	#extract T0 timepoints to calculate p-value differences
	if (skipT0 == FALSE){
		var0 = unlist(strsplit(dat,"\\."))
		var1= paste(var0[1:length(var0)-1],collapse=".")
		var2=paste(var1,"T0", sep='.') #baseline

		output_val[4] = t.test(group2[[dat]], group2[[var2]], paired = T, na.rm = T)$p.value
		output_val[5] = t.test(group3[[dat]], group3[[var2]], paired = T, na.rm = T)$p.value
		output_val[6] = t.test(all[[dat]], all[[var2]], paired = T, na.rm = T)$p.value
		
		output_val[7] = calc_effect_size(group2[[dat]], group2[[var2]])
		output_val[8] = calc_effect_size(group3[[dat]], group3[[var2]])
		output_val[9] = calc_effect_size(all[[dat]], all[[var2]])
		} else {
		output_val[4] = t.test(group2[[dat]], group1[[dat]], na.rm = T)$p.value
		output_val[5] = t.test(group3[[dat]], group1[[dat]], na.rm = T)$p.value
		output_val[6] = t.test(all[[dat]], group1[[dat]], na.rm = T)$p.value
		
		output_val[7] = calc_effect_size(group2[[dat]], group1[[dat]])
		output_val[8] = calc_effect_size(group3[[dat]], group1[[dat]])
		output_val[9] = calc_effect_size(all[[dat]], group1[[dat]])
		}
	
	if (sd(group2[[dat]], na.rm = TRUE) > 0 && sd(group3[[dat]], na.rm = TRUE) > 0){	
		#shapiro test to determine t test or wilcox
		n1 = shapiro.test(group2[[dat]])$p.value
		n2 = shapiro.test(group3[[dat]])$p.value
		
		#calculate magnitude
		mean_group2 = mean(group2[[dat]], na.rm=TRUE)
		mean_group3 = mean(group3[[dat]], na.rm=TRUE)
		magnitude = mean_group2/mean_group3
		if (magnitude < 1){			#convert to fold-change
			magnitude = 1/magnitude
			}
		output_val[2] = magnitude
		
		#calcuate effect size
		effect_size = calc_effect_size(group2[[dat]], group3[[dat]])
		output_val[3] = effect_size
		}
		
	#group effect
	p_one_tail = 1.0
	if (n1 > 0.05 && n2 > 0.05){
		p_one_tail = t.test(group2[[dat]], group3[[dat]], na.rm = T,alternative="two.sided")$p.value
		} else {		
		p_one_tail = wilcox.test(group2[[dat]], group3[[dat]], na.rm = T,alternative="two.sided")$p.value
		}
	p = round(p_one_tail, 4)
	output_val[1] = p
	
	if(draw_fig == TRUE){
		#figure parameters
		png(paste(pathtag,"fig_", filetag,"/","fig_",dat,".png", sep=''), width = 5, height = 3.3, units = 'in', res = 300)
		par(mfrow=c(1,1))
		par(ps = 9, cex = 1.0, cex.main = 1.0, cex.lab = 1.0, cex.axis = 1.0, cex.sub = 1.0)
		par(mai = c(0.8, 0.8, 0.3, 0.3))
		par(mgp = c(2.8,0.7,0))
		par(las=1)
		ylabel = dat
		
		ymin = max(0, min(plot_data$dat, na.rm = TRUE) - IQR(plot_data$dat, na.rm = T))
		ymax = median(plot_data$dat, na.rm = TRUE) + 2.5*IQR(plot_data$dat, na.rm = T)
		 # ymin = 0
		 # ymax = 1.5*IQR(plot_data$dat, na.rm = T)

		#boxplot(dat ~ Cohort, data = plot_data, at=c(2,3,1), outline = FALSE, ylim = c(ymin, ymax), ylab = ylabel, xlab = '', na.rm = T)
		#beeswarm(dat ~ Cohort, data = plot_data, at=c(2,3,1),pch = 16, ylim = c(ymin, ymax), ylab = ylabel, xlab = '', add = TRUE, pwcol = dot_col[as.numeric(Dose)], na.rm = T)
		#legend("topleft", legend = c("ALFA","ALFQ","ALFQA"), pch = 16, col = dot_col[1:3])
		
		boxplot(dat ~ Adjuvant, data = plot_data, at=c(2,3,4,1), outline = FALSE, ylim = c(ymin, ymax), ylab = ylabel, xlab = '', na.rm = T)
		beeswarm(dat ~ Adjuvant, data = plot_data, at=c(2,3,4,1),pch = 16, ylim = c(ymin, ymax), ylab = ylabel, xlab = '', add = TRUE, pwcol = dot_col[as.numeric(Dose)], na.rm = T)
		legend("topleft", legend = c("low dose", "high dose", "None"), pch = 16, col = dot_col[1:3])

		if (is.nan(output_val[1]) == F){
			#text(2, 0.9*ymax, paste("dose-effect p = ", round(output_val[1],3), sep=' '))
			text(2.5, 0.9*ymax, paste("adj-effect p = ", round(output_val[1],3), sep=' '))
			}
			
		min_p = min(output_val[4:6])
		if (is.nan(min_p) == F){
			text(2.5, 1*ymax, paste("vaccine-induced p = ", round(min_p,3), sep=' '))
			}
		

		dev.off()
		}
		
	output_val
}

run_analysis = function(filetag, pathtag, draw_fig, skipT0, assayID, exclude_col){
	info = read.table("./data/info.csv", header = T, sep =',')
	results = read.table(paste("./data/", filetag, ".csv", sep=''), header = T, sep =',')
	
	d1 = merge(info, results, by="Subject_ID")
	d1 = add_assayID(d1, assayID, exclude_col)	
	
	cols = colnames(d1)
	assayID = array(0.0, c(length(cols)))
	timepoint = array(0.0, c(length(cols)))
	group.p = array(1.0, c(length(cols)))
	group.magnitude = array(0.0, c(length(cols)))
	group.effect_size = array(0.0, c(length(cols)))
	vacc.grp2.p = array(1.0, c(length(cols)))
	vacc.grp3.p = array(1.0, c(length(cols)))
	vacc.all.p = array(1.0, c(length(cols)))
	vacc.grp2.effect_size = array(1.0, c(length(cols)))
	vacc.grp3.effect_size = array(1.0, c(length(cols)))
	vacc.all.effect_size = array(1.0, c(length(cols)))
	 
	for(i in 1:length(cols)){
		skip = F
		for(j in 1:length(exclude_col)){
			if (cols[i] == exclude_col[j]){
				skip = T
				}
			}
		if (skip == F){
			calc_val = plot_chart(cols[i], d1, filetag, pathtag, draw_fig, skipT0)
			group.p[i] = calc_val[1]
			group.magnitude[i] = calc_val[2]
			group.effect_size[i] = calc_val[3]
			vacc.grp2.p[i] = calc_val[4]
			vacc.grp3.p[i] = calc_val[5]
			vacc.all.p[i] = calc_val[6]
			vacc.grp2.effect_size[i] = calc_val[7]
			vacc.grp3.effect_size[i] = calc_val[8]
			vacc.all.effect_size[i] = calc_val[9]
			timepoint[i] = extract_timepoint(cols[i])
			assayID = extract_assayID(cols[i])
			}
		}

	#output_table = data.frame(cols, timepoint, pval, magnitude, effect_size, p_grp2, p_grp3, p_all, es_grp2, es_grp3, es_all)
	output_table = data.frame(cols, assayID, timepoint, 
								vacc.grp2.p, vacc.grp2.effect_size, 
								vacc.grp3.p, vacc.grp3.effect_size,
								vacc.all.p, vacc.all.effect_size,
								group.p, group.magnitude, group.effect_size)
								
	for(i in 1:length(exclude_col)){
		output_table <- output_table[!(output_table$cols==exclude_col[i]),]
		}
	#write.table(output_table, paste(pathtag,filetag, ".pvalue.csv", sep=''), sep=';', row.names = F)	
	output_table
	}

calc_qvalue = function(d1, filetag, pathtag, exclude_col){
	#d1 = read.table(paste(d1, pathtag,filetag, ".pvalue.csv", sep=''), header = T, sep =';')
	
	new_d1 = d1
	for(i in 1:length(exclude_col)){
		for(j in 1:length(d1$cols)){
			if (d1$cols[j] == exclude_col[i]){
				new_d1 = d1[-j,]
				}
			}
		d1 = new_d1
		}
	time = unique(d1$timepoint)
	output_table = "NA"

	for (i in 1:length(time)){
		sub_d1=subset(d1, timepoint == time[i]) #calculates q-values within each time point
		n_param = length(sub_d1$group.p)
		
		#calculate q-values for group-level differences 
		sub_d1$group.q = p.adjust(sub_d1$group.p, "BH", n_param)
		
		#calculate q-values for vaccine-induced effect (vs T0)
		sub_d1$vacc.grp2.q = array('NA', dim=c(n_param))
		sub_d1$vacc.grp3.q = array('NA', dim=c(n_param))
		sub_d1$vacc.all.q = array('NA', dim=c(n_param))
		
		if(is.na(sub_d1$vacc.all.p[1]) == FALSE){#for some measures, this is not calculable
			sub_d1$vacc.grp2.q = p.adjust(sub_d1$vacc.grp2.p, "BH", n_param)
			sub_d1$vacc.grp3.q = p.adjust(sub_d1$vacc.grp3.p, "BH", n_param)
			sub_d1$vacc.all.q= p.adjust(sub_d1$vacc.all.p, "BH", n_param)
			}
			
		if (i == 1){
			output_table = sub_d1
			} else {
			output_table = rbind(output_table, sub_d1)
			}
		}
	write.table(output_table, paste(pathtag,filetag, ".qvalue.csv", sep=''), sep=';', row.names = F, quote=FALSE)	
	}	
	
#DESCRIPTIVE STATISTICS FUNCTIONS END

#PCA ANALYSIS FUNCTIONS
plot_pca = function(pathtag, filename, imm_val, imm_cohort, dot_col){

	ir.cohort = imm_cohort
	
	#fill in missing values
	for (i in 1:length(imm_val)){	#columns
		for (j in 1:length(imm_val[,1])){	#subjects
			if (is.na(imm_val[j,i])){
				imm_val[j,i] = median(imm_val[,i], na.rm = TRUE)
				}
			}
		}
	
	ir.pca = prcomp(imm_val, center = TRUE, scale = TRUE)
	#print(summary(ir.pca))
	
	png(paste(pathtag, filename, sep=''), width = 5, height = 5, units = 'in', res = 300)
	
	#biplot(ir.pca, scale = 0.5)
	#plot(ir.pca, type ='l')
	
	g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.cohort, ellipse = TRUE, 
              circle = TRUE)
	g <- g + scale_color_discrete(name = '')
	g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
	dev.off()
}
#PCA ANALYSIS FUNCTIONS END

