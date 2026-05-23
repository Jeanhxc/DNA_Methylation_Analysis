#CG
callusCG= read.table("common_region_CG.txt",sep="\t",header=T)
cgtrun_selected_columns <- c("Chr","start","end","Transformed_1","Transformed_2","Transformed_3","Untransformed_1","Untransformed_2","Untransformed_3")
CGtrun <- callusCG[, cgtrun_selected_columns]

cutoff_CG=0.25

CGtrun[[4]] = CGtrun[[4]] + 0.0001
CGtrun[[7]] = CGtrun[[7]] + 0.0001

CGtrun$p_value=apply(CGtrun[,4:8],1,function(x){t.test(x[1:3], x[4:6], alternative="two.sided")$p.value})
CGtrun[[4]]=CGtrun[[4]]-0.0001
CGtrun[[7]]=CGtrun[[7]]-0.0001

CGtrun$Un_mean=apply(CGtrun[,7:9],1, mean)
CGtrun$Tr_mean=apply(CGtrun[,4:6],1, mean)
CGtrun$delta=CGtrun$Tr_mean-CGtrun$Un_mean

CG_trun=CGtrun[abs(CGtrun$delta)>=cutoff_CG&(CGtrun$p_value< 0.05),]
CG_trun$type <- ifelse (CG_trun$delta>0,"hyper","hypo")


