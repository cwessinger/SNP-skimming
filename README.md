# SNP-skimming
Python scripts to perform SNP-skimming

## Important notes

- These scripts require Python 2, NumPy, SciPy, and Pandas.

- The scripts assume a single population appears in the VCF file and all individual samples in the VCF are being analyzed. In this case, one provides the number of individuals as an argument to the script. 

  If a user has a more complicated scenario where the individuals to be analyzed are a subset of individuals that appear in the VCF file, one should modify the `indivRange` variable within the script to indicate the columns in which the samples to be analyzed appear in the VCF file.

- The scripts take as input a VCF file that has been filtered to the user's specification, in terms of allowed amount of missing data or range of allele frequency.  
  - For help on generating a VCF file, see **Preparing your VCF file** below.  
  - For help filtering a VCF file based on mapping quality, missing data, and very roughly estimated allele frequency, see **Filtering VCF file** below.


## Step 1: Estimate allele frequencies using MLEq.n.filter.py

Use **MLEq.n.filter.py** to find maximum likelihood estimator for frequency of the ref allele (_q_) for each site in the VCF file.  

In addition, this script discards sites that are actually fixed for the ref or the alt allele. Sometimes these sites sneak through despite filtering. (Basically, if the maximum likelihood estimate for *q* is 0 or 1, the site is discarded).

#### Usage

> python MLEq.n.filter.py [input VCF] [output prefix] [# indivs]

Arguments:
* input VCF: this is the VCF you want to use for your analyses
* output prefix: prefix for output filenames
* # indivs: number of individuals with data in the VCF file

#### Example

> `python MLEq.n.filter.py snps1.1.vcf snps1.2 291`

This creates four output files: 
1. `<prefix>.poly.vcf` <- A new .vcf file containing lines for non-fixed sites. (Header lines are removed)
1. `<prefix>.fixed.vcf` <- A new .vcf file containing lines for fixed sites. (Header lines removed).
1. `<prefix>.poly.Qs.txt` <- A tab-delimited file containing allele frequency estimates for ref base for non-fixed SNPs.
1. `<prefix>.fixed.Qs.txt` <- A tab-delimited file containing allele frequency estimates for ref base for fixed SNPs.

The `<prefix>.poly.vcf` and `<prefix>.poly.Qs.txt` files are used in downstream analyses, but users may also find the other two files useful to have.

## Step 2: Estimate τk values across all sites and all individuals using MLEtk.allindivs.py

Use **MLEtk.allindivs.py** to find maximum likelihood estimator for τk, the probability of detecting a heterozygote for read depth k. This estimate is based on data across all sites and all individuals for a given read depth. By default, the script produces a file of estimated τk values for k ranging from 1-15 reads. However, one could extend this to higher read depths if τk values in this range have not yet attained 1.

The input files should be `<prefix>.poly.vcf` and `<prefix>.poly.Qs.txt` files output by **MLEq.n.filter.py**.

#### Usage

> python MLEtk.allindivs.py [input VCF] [input qs file] [output τks file] [# indivs]

#### Example

> `python MLEtk.allindivs.py snps1.2.poly.vcf snps1.2.poly.Qs.txt snps1.2.poly.tks.txt 291`

The output file contains the number of sites, proportion heterozygous genotypes, and maximum likelihood estimate of τk (for k values ranging from 1–15). Depending on how this these results look, one could extend the calculations to higher read depth values by modifying the `max_kvalue` variable in the script.

## Step 2b: Estimate individual-specific τk values using MLEtk.byindivs.py

Before proceeding to step 3, it is advisable to examine how τk values vary across individuals. One can proceed with a global set of τk values (across all individuals - estimated using MLEtk.allindivs.py) only if τk values are fairly similar across individuals. 

Use **MLEtk.byindivs.py** to find τk values for each individual.

#### Usage

> python MLEtk.byindivs.py [input VCF] [input qs file] [output τks file] [# indivs]

#### Example

> `python MLEtk.byindivs.py snps1.2.poly.vcf snps1.2.poly.Qs.txt snps1.2.poly.tks.byindivs.txt 291`

## Step 3: Re-estimate allele frequencies (*q*s) using MLEq.reest.py.

Use **MLEq.reest.py** to re-estimate frequency of the reference allele, given read depth and global τk values estimated in step 2. 

#### Usage

> python MLEq.reest.py [input VCF] [input qs file] [input τks file] [output re-est qs file] [# indivs]

#### Example

> `python MLEq.reest.py snps1.2.poly.vcf snps1.2.poly.Qs.txt snps1.2.poly.tks.txt snps1.2.poly.reestQs.txt 291`

## Step 4: Estimate associations with phenotype using MLE 

Use **MLE.phenoassoc.py** to estimate associations between sites in the VCF file and a focal phenotype. 

Required input files:
1. VCF file for variable sites (output by step 1)
1. File containing re-estimated allele frequencies (output by step 3). Each line corresponds to a site in the VCF file.
1. Estimated global τk values (output by step 2)
1. Tab delimited file containing phenotypic values and a header line containing phenotype names. This file is read in by pandas read_table(). 

This script optimizes multiple parameters and can be time-consuming for large VCF files. I recommend splitting the VCF file (and corresponding allele frequency file) into many subsamples and running parallel jobs, then concatenating the output files for analysis.

#### Usage

> python MLE.phenoassoc.py [input VCF] [input re-est qs file] [input τks file] [input phenos file] [# indivs] [pheno name] [VCF ID]

Here, the [pheno name] argument should be the column header in the phenotype file for the focal phenotype to be analyzed.
The [VCF ID] argument should be the identifier (e.g., an integer) for subsample VCF file. This VCF ID is incorporated into the output file name. If one prefers not to split the file into parallel jobs, just put a placeholder here.

> `python MLE.phenoassoc.py snps1.2.poly.sub.0.vcf snps1.2.reestQs.sub.0.txt snps1.2.poly.tks.txt phenos.txt 291 stamenL 0`

This file outputs the following data for each SNP in the VCF (MLE = maximum likelihood estimate):
* `contig`: identifier of contig or scaffold in ref genome
* `position`: position of SNP
* `count`: number of individuals with data for the SNP
* `q`: allele frequency (read in from the input re-est qs file)
* `h0_mu`: MLE for sample mean phenotype under null hypothesis that phenotype does not depend on genotype
* `h0_sigma`: MLE for sample variance in phenotype under null hypothesis
* `h1_muRR`: MLE for mean phenotype for ref allele homozygotes under alternative hypothesis that genotypes differ in mean phenotype
* `h1_a`: MLE for additive effect of alt allele under alternative hypothesis
* `h1_sigma`: MLE for variance in phenotype for a given genotype under alternative hypothesis
* `LRT`: likelihood ratio test statistic for comparing the alternative hypothesis to the null hypothesis 
* `pval`: significance of the likelihood ratio test statistic (using chi-square distribution with 1 df)

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
