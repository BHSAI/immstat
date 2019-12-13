#filetag = c("ICC", "IGG", "AELI", "BELI", "BCF", "CYT")
filetag = c("ICC", "IGG", "BELI", "BCF")
#filetag = c("SERO")
pathtag = "./results_group_qval/" ##FILEPATH FOR STATUS OR PROTECTION
tag = ".group"

clust_list = read.table("ClustID.all.txt", header = T, sep=';')

for (i in 1:length(filetag)){
	print(paste(pathtag,filetag[i], ".qvalue.csv", sep=''))
	d1=read.table(paste(pathtag,filetag[i], ".qvalue.csv", sep=''), header = T, sep=';')
	time_list = unique(d1$timepoint)
	
	#initialize array
	tmp1 = subset(d1, timepoint == time_list[1])
	cols = tmp1$cols
	
	output_data = data.frame(cols)
	
	tp = c("T2", "T4", "T5", "T6", "T7")
	for(j in 1:length(tp)){
		if(is.na(match(tp[j], time_list))==FALSE){
			tmp = subset(d1, timepoint == tp[j])
			clustID = paste(tp[j],".clustID",sep='')
			pname = paste(tp[j],tag,".p",sep='')
			qname = paste(tp[j],tag,".q",sep='')
			vpname = paste(tp[j],".vacc.p",sep='')
			vqname = paste(tp[j],".vacc.q",sep='')
			v_class = paste(tp[j],".vacc.class", sep='') #vaccine-induced response classification
			g_class = paste(tp[j], tag,".class", sep='') #group-level differences classification
	
			min_vacc_p = array(1.0, dim=c(length(tmp$cols)))
			min_vacc_q = array(1.0, dim=c(length(tmp$cols)))
			vacc_class = array(0, dim=c(length(tmp$cols)))
			group_class = array(0, dim=c(length(tmp$cols)))
			clusterID = array(999, dim=c(length(tmp$cols))) #only vaccine-induced measures are clustered
			
			for(k in 1:length(tmp$cols)){
				min_vacc_p[k] = min(c(tmp$p_grp2[k], tmp$p_grp3[k], tmp$p_all[k]))
				min_vacc_q[k] = min(c(tmp$q_grp2[k], tmp$q_grp3[k], tmp$q_all[k]))
				if(min_vacc_p[k] < 0.05 & min_vacc_q[k] < 0.20){
					vacc_class[k] = 1
					if(tmp$pval[k] < 0.05 && tmp$q_val[k] < 0.20){
						group_class[k] = 1
						}
					}
					
				time_col = paste(tmp$cols[k],tp[j],sep='.')
				print(time_col)
				find_col = match(time_col, clust_list$cols)
				if(is.na(find_col)==FALSE){
					clusterID[k] = clust_list$clust_ID[find_col]
					}
				}
			
			output_data[[vpname]] = min_vacc_p
			output_data[[vqname]] = min_vacc_q
			output_data[[pname]] = tmp$pval
			output_data[[qname]] = tmp$q_val	

			output_data[[v_class]] = vacc_class
			output_data[[g_class]] = group_class
			output_data[[clustID]] = clusterID
			}
		}
		
	write.table(output_data, paste(pathtag,filetag[i], ".merged.csv", sep=''), sep=';', row.names = F)
	}
	
	
	
	
	
	

