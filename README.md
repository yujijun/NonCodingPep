# NonCodingPep
This is a neoantigen prediction pipeline which include peptides from non-coding region 

# Cancer vs normal pipeline code
# Step 1
../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode genomeGenerate      --genomeDir C_genome_dir      --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa      --runThreadN 6
../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode genomeGenerate      --genomeDir N_genome_dir      --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa      --runThreadN 6

# Step 2
 nohup ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeDir C_genome_dir       --readFilesIn C_1_paired.fq.gz C_2_paired.fq.gz       --runThreadN 8       --readFilesCommand zcat       --outFilterType BySJout       --outFileNamePrefix C_ &
nohup ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeDir N_genome_dir       --readFilesIn N_1_paired.fq.gz N_2_paired.fq.gz       --runThreadN 8       --readFilesCommand zcat       --outFilterType BySJout       --outFileNamePrefix N_ &

# Step 3
    nohup  ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
     --genomeDir C_genome_dir \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbFileChrStartEnd C_SJ.out.tab \
     --sjdbOverhang 75 \
     --runThreadN 8 &

    nohup  ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
     --genomeDir N_genome_dir \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbFileChrStartEnd N_SJ.out.tab \
     --sjdbOverhang 75 \
     --runThreadN 8 &

# Step 4
   nohup   ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeDir C_genome_dir \
     --readFilesIn C_1_paired.fq.gz C_2_paired.fq.gz \
     --runThreadN 8 \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:ENCSR000COQ1 LB:library PL:illumina PU:machine SM:GM12878 \
     --outFileNamePrefix C_ &
   
  ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeDir N_genome_dir \
     --readFilesIn N_1_paired.fq.gz N_2_paired.fq.gz \
     --runThreadN 8 \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:ENCSR000COQ1 LB:library PL:illumina PU:machine SM:GM12878 \
     --outFileNamePrefix N_

# Step 5.   
samtools index C_Aligned.sortedByCoord.out.bam
samtools index N_Aligned.sortedByCoord.out.bam

# Step 6. 
mkdir tmp 
java -Xmx120G -Djava.io.tmpdir=tmp/ -jar ../../software/GenomeAnalysisTK.jar -T SplitNCigarReads \
                  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -I C_Aligned.sortedByCoord.out.bam \
                  -o C_split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS 
(samtools view -H C_split.bam; samtools view C_split.bam| grep -w 'NH:i:1') \
  | samtools view -Sb -  > C_final.uniq.bam

mkdir tmp 
java -Xmx120G -Djava.io.tmpdir=tmp/ -jar ../../software/GenomeAnalysisTK.jar -T SplitNCigarReads \
                  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -I N_Aligned.sortedByCoord.out.bam \
                  -o N_split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS 
(samtools view -H N_split.bam; samtools view N_split.bam| grep -w 'NH:i:1') \
  | samtools view -Sb -  > N_final.uniq.bam

# Step 7. 
(samtools view -H N_split.bam; samtools view N_split.bam| grep -w 'NH:i:1') \
  | samtools view -Sb -  > N_final.uniq.bam

(samtools view -H C_split.bam; samtools view C_split.bam| grep -w 'NH:i:1') \
  | samtools view -Sb -  > C_final.uniq.bam

# Step 8. 
samtools index C_final.uniq.bam
samtools index N_final.uniq.bam

# Step 9. 
normal1_pileup='samtools mpileup -q 1 -f Homo_sapiens.GRCh38.dna.primary_assembly.fa C_final.uniq.bam';
tumor1_pileup='samtools mpileup -q 1 -f Homo_sapiens.GRCh38.dna.primary_assembly.fa N_final.uniq.bam';
sample=test
java -Xmx120G -Djava.io.tmpdir=tmp/ -jar /mnt/data/meng/software/VarScan.v2.4.4.jar processSomatic test.snp

# Cancer vs reference pipeline code
1. data download
# sample：
normal: SRR8668622/SRR8668621/SRR8668620/SRR8668619
cancer: SRR7094748

# reference genome: 
nohup wget ftp://ftp.ensembl.org/pub/release102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.chr.gtf.gz & (人类参考基因组注释文件)
nohup wget ftp://ftp.ensembl.org/pub/release102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz & (人类参考基因组)
nohuo wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa &
（#组装版本，The STAR manual tells us that Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz is an acceptable file to use）

2. Data preprocessing -- from sra to fq

# 软件下载安装  sratoolkit.2.10.9
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.2.10.9-ubuntu64.tar.gz

# 确认下载数据为单端还是双端（返回值为4，为单端SE；若返回值为8，则为双端PE）
../../software/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump -X 1 --split-spot -Z SRR7094748.1 | wc -l

# 格式转换
for i in SRR*; do echo $i; ../../sratoolkit.2.10.9-ubuntu64/bin/fastq-dump 

3. Data preprocessing -- quality control

# 软件下载安装  FastQC v0.11.9 (Linux zip file)
              unzip fastqc_v0.11.9.zip
              chmod 755 fastqc
# 质控
              ../../software/FastQC/fastqc C_1.fastq
              ../../software/FastQC/fastqc C_2.fastq
              ../../software/FastQC/fastqc N_1.fastq
              ../../software/FastQC/fastqc N_2.fastq

4. Data preprocessing -- filtering

#  软件下载安装    
unzip Trimmomatic-0.39.zip

#  过滤脚本
* for Cancer sample
java -jar /mnt/data/meng/software/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE \
C_1.fastq C_2.fastq \
C_1_paired.fq.gz C_1_unpaired.fq.gz C_2_paired.fq.gz C_2_unpaired.fq.gz \
ILLUMINACLIP://mnt/data/meng/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

* for Normal sample
java -jar /mnt/data/meng/software/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE \
N_1.fastq N_2.fastq \
N_1_paired.fq.gz N_1_unpaired.fq.gz N_2_paired.fq.gz N_2_unpaired.fq.gz \
ILLUMINACLIP://mnt/data/meng/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

5. Data preprocessing -- re quality control

# 再次质控
          ../../software/FastQC/fastqc C_1_paired.fq.gz
          ../../software/FastQC/fastqc C_2_paired.fq.gz
          ../../software/FastQC/fastqc N_1_paired.fq.gz
          ../../software/FastQC/fastqc N_2_paired.fq.gz

6.  Data processing -- mapping

6.1 为GATK需求，首先产生基因组索引文件

# samtools fadix 
        samtools-1.11.tar.bz2
        tar jxvf samtools-1.9.tar.bz2
        cd samtools-1.11
        ./configure --prefix=/where/to/install
        make
        make install
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

#  picard
        https://github.com/broadinstitute/picard/releases/tag/2.25.0
        java -jar ../../software/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.primary_assembly.fa O=Homo_sapiens.GRCh38.dna.primary_assembly.dict

#  mapping -- STAR
    
     STAR （https://github.com/alexdobin/STAR/releases/tag/2.7.7a）
     STAR-2.7.7a.tar.gz
     tar -zxvf STAR-2.7.7a.tar.gz
    
     ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
     --genomeDir C_genome_dir \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --runThreadN 8 

     * 1st pass
     ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeDir C_genome_dir \
      --readFilesIn C_1_paired.fq.gz C_2_paired.fq.gz \
      --runThreadN 8 \
      --readFilesCommand zcat \
      --outFilterType BySJout \
      --outFileNamePrefix C_  
    
    * index
     ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
     --genomeDir C_genome_dir \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbFileChrStartEnd C_SJ.out.tab \
     --sjdbOverhang 75 \
     --runThreadN 8 
     
     2st pass
     ../../software/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeDir C_genome_dir \
     --readFilesIn C_1_paired.fq.gz C_2_paired.fq.gz \
     --runThreadN 8 \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:ENCSR000COQ1 LB:library PL:illumina PU:machine SM:GM12878 \
     --outFileNamePrefix C_
     (The STAR manual tells us that we need make 2-pass mapping)

     samtools index C_Aligned.sortedByCoord.out.bam

7.  Data processing -- variants calling
         
         https://storage.cloud.google.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
         mkdir tmp 
         java -Xmx120G -Djava.io.tmpdir=tmp/ -jar ../software/GenomeAnalysisTK.jar -T SplitNCigarReads \
                  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -I C_Aligned.sortedByCoord.out.bam \
                  -o C_split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS 

         (samtools view -H C_split.bam; samtools view C_split.bam| grep -w 'NH:i:1') \
         | samtools view -Sb -  > C_final.uniq.bam

         samtools index C_final.uniq.bam

         ls C_final.uniq.bam  > C_bam.list
         java -Xmx120G -Djava.io.tmpdir=tmp/ -jar ../../software/GenomeAnalysisTK.jar -T HaplotypeCaller \
                  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa  -I C_bam.list \
                  -dontUseSoftClippedBases \
                  -stand_call_conf 20.0 \
                  -o C_output.gatk.vcf.gz

         java -Xmx120G -Djava.io.tmpdir=tmp/ -jar ../software/GenomeAnalysisTK.jar -T VariantFiltration \
                  -R Homo_sapiens.GRCh38.dna.primary_assembly.fa -V C_output.gatk.vcf.gz \
                  -window 35 -cluster 3 \
                  -filterName FS -filter "FS > 30.0" \
                  -filterName QD -filter "QD < 2.0" \
                  -o C_final.vcf

8.  Annoation
         https://annovar.openbioinformatics.org/en/latest/user-guide/download/
         tar xvfz annovar.latest.tar.gz
         nohup ../../software/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar ensGene humandb/ &
         
         perl /data/software/annovar/table_annovar.pl C_final.vcf  \
         --outfile C_annovar \
         --buildver hg38 \
         --protocol refGene \
         --operation g \
         --vcfinput \
         humandb/