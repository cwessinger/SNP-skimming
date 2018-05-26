# SNP-skimming
Python scripts to perform SNP-skimming

## Important details regarding these python scripts

**Use python 2.7**

## Preparing your VCF file

SNP-skimming takes a .vcf file as input. 

Here we describe the steps we like to take to generate this file, but other methods could be used.  
(In this case, `refcontigs.fa` refers to reference genome assembly and `ind1.fq` is sequence data for one sample)  

1. **Start with a single .fastq file for each individual** that has been filtered and trimmed for quality and to remove adapter sequence.  One option is steps 1 & 2 of ipyrad -> http://ipyrad.readthedocs.io

1. **Map .fastq files to the reference genome assembly using BWA** -> https://sourceforge.net/projects/bio-bwa  
    `bwa mem refcontigs.fa trimmed.ind1.fq.gz > ind1.sam`

1. **Convert .sam to .bam using samtools** -> http://www.htslib.org  
    `samtools view -bS ind1.sam > ind1.unsorted.bam`

1. **Sort .bam using samtools**  
    `samtools sort ind1.unsorted.bam > ind1.bam`

1. **Add read group tags using picard** -> https://github.com/broadinstitute/picard  
    `java -jar picard.jar AddOrReplaceReadGroups I=ind1.bam O=ind1.RG.bam SO=coordinate RGID=SeqRUN# RGLB=ind1.bam RGPL=illumina RGPU=ind1.bam RGSM=ind1.bam VALIDATION_STRINGENCY=LENIENT`

1. **Index .bam file using samtools**  
    `samtools index ind1.RG.bam ind1.RG.bai`

1. **Call variants using GATK** -> https://software.broadinstitute.org/gatk/  
    `java -jar GenomeAnalysisTK.jar -R refcontigs.fa -T UnifiedGenotyper -I ind1.RG.bam -I ind2.RG.bam -I ind3.RG.bam -ploidy 2 -rf BadCigar -o data.vcf`
