#1. Quality Control
fastqc sample.fastq

#2.Remove PCR duplicate
BSseeker2/rmDupPE.pl sample_R1.fastq sample_R2.fastq sample_R1_rmdup.fastq.gz sample_R2_rmdup.fastq.gz

#3. reverse—complement
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


#4.Trim adaptor
java -jar trimmomatic-0.39.jar SE -phred33 sample_R1_rmdup.fastq.gz sample_R1_trimming.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#5.Alignment 
BSseeker2/bs_seeker2-build.py -f MSUv7.fa --aligner=bowtie2 -d ./BS2_bt2_Index
BSseeker2/bs_seeker2-align.py -i sample_R1_trimming.fastq.gz -g MSUv7.fa --aligner=bowtie2 -o sample_R1.bam -m 3 -d ./BS2_bt2_Index
samtools merge sample.bam sample_R1.bam sample_R2.bam

#6.Methylatino calling
BSseeker2/bs_seeker2-call_methylation.py -i sample.bam -o sample_out -d ~/tools/BSseeker2/BS2_bt2_Index/MSUv7.fa_bowtie2
