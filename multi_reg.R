d1 = read.table("mal71_summary.csv", header = T, sep =';')
#d1 = read.table("depletion_pool.csv", header = T, sep =';')
d1_sub = d1[,5:26]
#d1_sub = d1_sub * scale(d1_sub)

nsub = length(d1$Group)
dummy = array(0, c(nsub))

for(i in 1:nsub){
	if (d1$Group[i] == "RRRx"){
		dummy[i] = 1
		}
	}

d1_sub$dummy = dummy

fit_MFI_CS = lm(MFI2 ~ CS_Ig1 + CS_Ig2 + CS_Ig3 + CS_Ig4, data = d1_sub) 
fit_MFI_NANP = lm(MFI2 ~ NANP_Ig1 + NANP_Ig2 + NANP_Ig3 + NANP_Ig4, data = d1_sub) 
fit_MFI_PF16 = lm(MFI2 ~ PF16_Ig1 + PF16_Ig2 + PF16_Ig3 + PF16_Ig4, data = d1_sub) 

fit_MFI_CS_dummy = lm(MFI2 ~ CS_Ig1 + CS_Ig2 + CS_Ig3 + CS_Ig4 + dummy, data = d1_sub) 
fit_MFI_NANP_dummy = lm(MFI2 ~ NANP_Ig1 + NANP_Ig2 + NANP_Ig3 + NANP_Ig4 + dummy, data = d1_sub) 
fit_MFI_PF16_dummy = lm(MFI2 ~ PF16_Ig1 + PF16_Ig2 + PF16_Ig3 + PF16_Ig4 + dummy, data = d1_sub) 

fit_M2_CS = lm(M2 ~ CS_Ig1 + CS_Ig2 + CS_Ig3 + CS_Ig4, data = d1_sub) 
fit_M2_NANP = lm(M2 ~ NANP_Ig1 + NANP_Ig2 + NANP_Ig3 + NANP_Ig4, data = d1_sub) 
fit_M2_PF16 = lm(M2 ~ PF16_Ig1 + PF16_Ig2 + PF16_Ig3 + PF16_Ig4, data = d1_sub) 

fit_M2_CS_dummy = lm(M2 ~ CS_Ig1 + CS_Ig2 + CS_Ig3 + CS_Ig4 + dummy, data = d1_sub) 
fit_M2_NANP_dummy = lm(M2 ~ NANP_Ig1 + NANP_Ig2 + NANP_Ig3 + NANP_Ig4 + dummy, data = d1_sub) 
fit_M2_PF16_dummy = lm(M2 ~ PF16_Ig1 + PF16_Ig2 + PF16_Ig3 + PF16_Ig4 + dummy, data = d1_sub) 
