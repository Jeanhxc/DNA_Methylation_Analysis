import pandas as pd
import gzip 
import os
sample_files = ["sample_CGmap/Transformed_1_out.CGmap.gz","sample_CGmap/Transformed_2_out.CGmap.gz","sample_CGmap/Transformed_3_out.CGmap.gz",
                "sample_CGmap/Untransformed_1_out.CGmap.gz","sample_CGmap/Untransformed_2_out.CGmap.gz","sample_CGmap/Untransformed_3_out.CGmap.gz",]
Contexts = ["CG", "CHG","CHH"]

output_directory = "200bp_output"
for sample_file in sample_files:
    for Context in Contexts:
        base_name = os.path.basename(sample_file)
        sample_name = base_name.split("_")[0]  
        output_file = os.path.join(output_directory, "{}_{}.csv".format(sample_name, Context))
        CGmap = pd.read_csv(sample_file, sep="\t", names=["Chr", "Nucleotide", "Position", "Context", "Dinucleotide_context", "Level", "#_of_C", "#_of_C+T"])
        CGmap = CGmap[CGmap["#_of_C+T"] >= 4]
        CGmap = CGmap[CGmap["Context"] != "--"].reset_index().drop(columns="index")
        CGmap["Context"] = CGmap["Context"].astype(str)
        CGmap = CGmap[CGmap["Context"].str.contains(Context)]
        CGmap = CGmap.drop(["Nucleotide", "Dinucleotide_context", "#_of_C", "#_of_C+T"], 1)

        CGmap["Bin"] = pd.cut(CGmap["Position"],[i*200 for i in range(max(CGmap["Position"])//200+2)],right=True)
        def extract_start(interval):
            return interval.left
        CGmap["start"] = CGmap["Bin"].apply(extract_start)
        CGmap["start"] = CGmap["start"].astype(str)
        CGmap["start"] = pd.to_numeric(CGmap["start"].str.replace('(', ''))
        CGmap["start"] = CGmap["start"].astype(int)

        CGmap["Count"]=CGmap.groupby(["Chr","Context","Bin"])["start"].transform("count")
        condition=CGmap["Count"] > 4 
        CGmap=CGmap[condition] 
        CGmap=CGmap.drop(["start","Count"],axis=1)
        CGmap = CGmap.groupby(["Chr","Context","Bin"],as_index=False)[["Level"]].mean()
        CGmap["Level"]=(CGmap["Level"]).round(6)
        CGmap['Level'] = CGmap['Level'].fillna("-")
    
        CGmap.to_csv(output_file, sep='\t', index=False, header=True)


#CG_common_region_file
Transformed_1=pd.read_csv("/200bp_output/Transformed_1_CG.csv",sep="\t", skiprows=1, na_values=["-"]).dropna()
Transformed_2=pd.read_csv("/200bp_output/Transformed_2_CG.csv",sep="\t", skiprows=1, na_values=["-"]).dropna()
Transformed_3=pd.read_csv("/200bp_output/Transformed_3_CG.csv",sep="\t", skiprows=1, na_values=["-"]).dropna()
Untransformed_1=pd.read_csv("/200bp_output/Untransformed_1_CG.csv",sep="\t", skiprows=1, na_values=["-"]).dropna()
Untransformed_2=pd.read_csv("/200bp_output/Untransformed_2_CG.csv",sep="\t", skiprows=1, na_values=["-"]).dropna()
Untransformed_3=pd.read_csv("/200bp_output/Untransformed_3_CG.csv",sep="\t", skiprows=1, na_values=["-"]).dropna()

commen=pd.merge(Transformed_1,Transformed_2,how="inner",on=['Chr','Context','Bin'])
commen=pd.merge(commen,Transformed_3,how="inner",on=['Chr','Context','Bin'])
commen=pd.merge(commen,Untransformed_1,how="inner",on=['Chr','Context','Bin'])
commen=pd.merge(commen,Untransformed_2,how="inner",on=['Chr','Context','Bin'])
commen=pd.merge(commen,Untransformed_3,how="inner",on=['Chr','Context','Bin'])

commen["start"]=commen["Bin"].map(lambda x:x.split(",")[0])
commen["end"]=commen["Bin"].map(lambda x:x.split(",")[1])
commen["start"] = pd.to_numeric(commen["start"].str.replace('(', ''))
commen["end"] = pd.to_numeric(commen["end"].str.replace(']', ''))
commen=commen[["Chr","start","end","Transformed_1","Transformed_2","Transformed_3","Untransformed_1","Untransformed_2","Untransformed_3"]]
commen.to_csv("common_region_CG.txt",sep='\t',index=False,header=True)
