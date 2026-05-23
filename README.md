# DNA_Methylation_Analysis
Precess Pipeline
------------------------------
**Step1: Quality control**

Tool: FastQC
```bash
#wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip unzip fastqc_v0.11.8.zip
fastqc sample.fastq
#View result using web browser
```
<img width="1055" height="713" alt="image" src="https://github.com/user-attachments/assets/c7f3d587-070b-4936-958e-8038fa1ec8af" />



**Step2: Remove PCR duplicate**

Tool: BS-Seeker2
```bash
#wget https://paoyang.ipmb.sinica.edu.tw/meth_course/bowtie2-2.3.5.1-linux-x86_64.zip unzip bowtie2-2.3.5.1-linux-x86_64.zip
BSseeker2/rmDupPE.pl sample_R1.fastq sample_R2.fastq sample_R1_rmdup.fastq.gz sample_R2_rmdup.fastq.gz
```

**Step3: reverse—complement**

```bash
import gzip
with gzip.open('sample_R2_rmdup.fastq.gz') as reader, gzip.open('sample_R2_rmdup_reverse.fastq.gz', 'w') as writer:
    for index, line in enumerate(reader):
        if index % 4 == 1:
                string=str(line)
                reversed_string=string[::-1]
                comp_reversed_string=reversed_string.replace("A","t").replace("T","a").replace("G","c").replace("C","g")
                com_rev_str=comp_reversed_string.upper()
                writer.write(com_rev_str+"\n")
        elif index %4==2:
                writer.write(line.strip("\n"))
        elif index %4==3:
                qual=str(line)
                reversed_qual=qual[::-1]
                writer.write(reversed_qual+"\n")
        elif index %4==0:
                writer.write(line.strip("\n"))
```

**Step4: Trim adaptor**

Tool: Trimmomatic
```bash
#wget [https://paoyang.ipmb.sinica.edu.tw/meth_course/bowtie2-2.3.5.1-linux-x86_64.zip unzip bowtie2-2.3.5.1-linux-x86_64.zip]
java -jar trimmomatic-0.39.jar SE -phred33 sample_R1_rmdup.fastq.gz sample_R1_trimming.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#View result using web browser
```
<img width="1044" height="708" alt="image" src="https://github.com/user-attachments/assets/28fb8dbe-56cc-4e38-bba1-424cdd4fa739" />



**Step5: Alignment & Methylatino calling**

Tool: BS-Seeker2

**BS-Seeker2:align**
```bash
#Download bowtie2
#wget https://paoyang.ipmb.sinica.edu.tw/meth_course/bowtie2-2.3.5.1-linux-x86_64.zipunzip bowtie2-2.3.5.1-linux-x86_64.zip
#Download BS-Seeker2
#git clone https://github.com/BSSeeker/BSseeker2.git
BSseeker2/bs_seeker2-build.py -f MSUv7.fa --aligner=bowtie2 -d ./BS2_bt2_Index

BSseeker2/bs_seeker2-align.py -i sample_R1_trimming.fastq.gz -g MSUv7.fa --aligner=bowtie2 -o sample_R1.bam -m 3 -d ./BS2_bt2_Index
samtools merge sample.bam sample_R1.bam sample_R2.bam
```

**BS-Seeker2:call methylation**
```bash
BSseeker2/bs_seeker2-call_methylation.py -i sample.bam -o sample_out -d ~/tools/BSseeker2/BS2_bt2_Index/MSUv7.fa_bowtie2
zless sample_out.CGmap.gz
```
<img width="1059" height="658" alt="image" src="https://github.com/user-attachments/assets/07cfceb0-f34b-4b05-ad3a-4077041aa3c5" />




Indentift DMR Pipeline
------------------------------

```bash
callusCG= read.table("TNG67_TrUn_callus_common_region_CG.txt",sep="\t",header=T)
cgtrun_selected_columns <- c("Chr","start","end","T17","T32","T65","C02","C04")
CGtrun <- callusCG[, cgtrun_selected_columns]

cutoff_CG=0.25

CGtrun[[4]] = CGtrun[[4]] + 0.0001
CGtrun[[7]] = CGtrun[[7]] + 0.0001

CGtrun$p_value=apply(CGtrun[,4:8],1,function(x){t.test(x[1:3], x[4:5], alternative="two.sided")$p.value})
CGtrun[[4]]=CGtrun[[4]]-0.0001
CGtrun[[7]]=CGtrun[[7]]-0.0001

CGtrun$Un_mean=apply(CGtrun[,7:8],1, mean)
CGtrun$Tr_mean=apply(CGtrun[,4:6],1, mean)
CGtrun$delta=CGtrun$Tr_mean-CGtrun$Un_mean

CG_trun=CGtrun[abs(CGtrun$delta)>=cutoff_CG&(CGtrun$p_value< 0.05),]
CG_trun$type <- ifelse (CG_trun$delta>0,"hyper","hypo")

write.table(CG_trun,"DMR_trun_0.25_CG.txt",sep="\t",quote=F,col.names=T,row.names=F)

```
