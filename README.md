# SNP-skimming
Python scripts to perform SNP-skimming

## Important notes

* These scripts require Python 2, NumPy, SciPy, and Pandas.

* The scripts assume a single population appears in the VCF file and all individual samples in the VCF are being analyzed. In this case, one provides the number of individuals as an argument to the script. 

  If a user has a more complicated scenario where the individuals to be analyzed are a subset of individuals included in the VCF file, one should modify the `indivRange` variable within the script to indicate the columns in which the samples to be analyzed appear in the VCF file.

* The scripts take as input a VCF file that has been filtered to the user's specification, in terms of allowed amount of missing data or range of allele frequency.  
  - For help on generating a VCF file, see **Preparing your VCF file** below.  
  - For help filtering a VCF file based on mapping quality, missing data, and very roughly estimated allele frequency, see **Filtering VCF file** below.


## Step 1: Estimate allele frequencies using MLEq.n.filter.py

Use **MLEq.n.filter.py** to find maximum likelihood estimator for frequency of the ref allele (_q_) for each site in the VCF file.  

In addition, this script discards sites that are actually fixed for the ref or the alt allele. Sometimes these sites sneak through despite filtering. (Basically, if the maximum likelihood estimate for *q* is 0 or 1, the site is discarded).

### Usage

> python MLEq.n.filter.py [input VCF] [output prefix] [n indivs]  

- input VCF: this is the VCF you want to use for your analyses  
- output prefix: prefix for output filenames  
- n indivs: number of individuals with data in the VCF file  

### Example

```
python MLEq.n.filter.py toy.vcf toy 291
```

This creates four output files: 
- **\<prefix\>.poly.vcf** <- A new .vcf file containing lines for non-fixed sites. (Header lines are removed)  
- **\<prefix\>.fixed.vcf** <- A new .vcf file containing lines for fixed sites. (Header lines removed)  
- **\<prefix\>.poly.Qs.txt** <- A tab-delimited file containing allele frequency estimates for ref base for non-fixed SNPs.  
- **\<prefix\>.fixed.Qs.txt** <- A tab-delimited file containing allele frequency estimates for ref base for fixed SNPs.  

The **\<prefix\>.poly.vcf** and **\<prefix\>.poly.Qs.txt** files are used in downstream analyses, but users may also find the other two files useful to have.
 
The allele frequency output files contain the following information:
- contig: identifier of contig or scaffold in ref genome  
- position: position of SNP  
- indivs: number of individuals with data for the SNP  
- MLE_q: maximum likelihood estimate for reference allele frequency (*q*)  


## Step 2: Estimate τk values across all sites and all individuals using MLEtk.allindivs.py

Use **MLEtk.allindivs.py** to find maximum likelihood estimator for τk, the probability of detecting a heterozygote for read depth *k*. This estimate is based on data across all sites and all individuals for a given read depth. By default, the script produces a file of estimated τk values for *k* ranging from 1-15 reads. However, one could extend this to higher read depths if τk values hasn't reached 1 by *k*=15.

### Usage

> python MLEtk.allindivs.py [input VCF] [input qs file] [output τks file] [n indivs]

* input VCF: this should be the `<prefix>.poly.vcf` file output by **MLEq.n.filter.py**
* input qs file: this should be the `<prefix>.poly.Qs.txt` file output by **MLEq.n.filter.py**
* output τks file: desired filename for output
* n indivs: number of individuals with data in the VCF file

### Example

```
python MLEtk.allindivs.py toy.poly.vcf toy.poly.Qs.txt toy.poly.tks.txt 291
```

The output file contains the following data:
* depth: read depth (*k*)
* totsites: number of sites with depth *k*
* propRA: proportion of genotypes at depth *k* that appear heterozygous
* MLE_tk: maximum likelihood estimate for τk

Depending on how this these results look, one could extend the calculations to higher read depth values by modifying the `max_kvalue` variable in the script.

## Step 2b: Estimate individual-specific τk values using MLEtk.byindivs.py

Before proceeding to step 3, it is advisable to examine how τk values vary across individuals. One can proceed with a global set of τk values (across all individuals - estimated using MLEtk.allindivs.py) only if τk values are fairly similar across individuals. 

Use **MLEtk.byindivs.py** to find τk values for each individual.

### Usage

> python MLEtk.byindivs.py [input VCF] [input qs file] [output τks file] [n indivs]

* input VCF: this should be the `<prefix>.poly.vcf` file output by **MLEq.n.filter.py**
* input qs file: this should be the `<prefix>.poly.Qs.txt` file output by **MLEq.n.filter.py**
* output τks file: desired filename for output
* n indivs: number of individuals with data in the VCF file


### Example

```
python MLEtk.byindivs.py toy.poly.vcf toy.poly.Qs.txt toy.poly.tks.byindivs.txt 291
```

The output file contains the following data:
* indiv: individual sample index
* depth: read depth (*k*)
* totsites: number of sites with depth *k*
* propRA: proportion of genotypes at depth *k* that appear heterozygous
* MLE_tk: maximum likelihood estimate for τk

## Step 3: Re-estimate allele frequencies (*q*s) using MLEq.reest.py.

Use **MLEq.reest.py** to re-estimate frequency of the reference allele, given read depth and global τk values estimated in step 2. 

### Usage

> python MLEq.reest.py [input VCF] [input qs file] [input τks file] [output re-est qs file] [n indivs]

* input VCF: this should be the `<prefix>.poly.vcf` file output by **MLEq.n.filter.py**
* input qs file: this should be the `<prefix>.poly.Qs.txt` file output by **MLEq.n.filter.py**
* input τks file: this is the τks file output by **MLEtk.allindivs.py**
* output re-est qs file: desired output filename
* n indivs: number of individuals with data in the VCF file


### Example
```
python MLEq.reest.py toy.poly.vcf toy.poly.Qs.txt toy.poly.tks.txt toy.reestQs.txt 291
```

The output file contains the following data:
* contig: identifier of contig or scaffold in ref genome
* position: position of SNP
* indivs: number of individuals with data for the SNP
* initial_q: initial maximum likelihood estimate for reference allele frequency (*q*) from input qs file
* reestimate_q: new ML estimate for *q*


## Step 4: Estimate associations with phenotype using MLE 

Use **MLE.phenoassoc.py** to estimate associations between sites in the VCF file and a focal phenotype. 

This script optimizes multiple parameters and can be time-consuming for large VCF files. **I recommend splitting the VCF file (and corresponding allele frequency file) into many subsamples and running parallel jobs,** then concatenating the output files for analysis.

### Usage

> python MLE.phenoassoc.py [input VCF] [input re-est qs file] [input τks file] [input phenos file] [n indivs] [pheno name] [VCF ID]

* input VCF: this should be the `<prefix>.poly.vcf` file output by **MLEq.n.filter.py**
* input re-est qs file: this should be the `<prefix>.poly.Qs.txt` file output by **MLEq.reest.py**
* input τks file: this is the τks file output by **MLEtk.allindivs.py**
* input phenos file: tab-delimited file containing phenotypic values and a header line containing phenotype names. THis file is read by Pandas read_table(). See file **phenos.txt** for an example file.
* output re-est qs file: desired output filename
* n indivs: number of individuals with data in the VCF file
* pheno name: name of phenotype to be name of the focal phenotype to be analyzed and should match the column header in the phenotype file.
* VCF ID: the identifier (e.g., an integer) for subsample VCF file. This VCF ID is incorporated into the output file name. If one prefers not to split the file into parallel jobs, just put a placeholder here.

### Example

```
python MLE.phenoassoc.py toy.poly.vcf toy.reestQs.txt toy.poly.tks.txt phenos.txt 291 StamenL 0
```

This file outputs the following data for each SNP in the VCF (MLE = maximum likelihood estimate):
* contig: identifier of contig or scaffold in ref genome
* position: position of SNP
* count: number of individuals with data for the SNP
* q: allele frequency (read in from the input re-est qs file)
* h0_mu: MLE for sample mean phenotype under null hypothesis that phenotype does not depend on genotype
* h0_sigma: MLE for sample variance in phenotype under null hypothesis
* h1_muRR: MLE for mean phenotype for ref allele homozygotes under alternative hypothesis that genotypes differ in mean phenotype
* h1_a: MLE for additive effect of alt allele under alternative hypothesis
* h1_sigma: MLE for variance in phenotype for a given genotype under alternative hypothesis
* LRT: likelihood ratio test statistic for comparing the alternative hypothesis to the null hypothesis 
* pval: significance of the likelihood ratio test statistic (using chi-square distribution with 1 df)

# Preparing your VCF file

SNP-skimming takes a .vcf file as input. 

Here I describe the steps we like to take to generate this file, but other methods could be used.  
(In this case, `refcontigs.fa` refers to reference genome assembly and `ind1.fq` is sequence data for one sample)  

1. **Start with a single .fastq file for each individual** that has been filtered and trimmed for quality and to remove adapter sequence.  One option is steps 1 & 2 of [ipyrad](http://ipyrad.readthedocs.io)

2. **Map .fastq files to the reference genome assembly using [BWA](https://sourceforge.net/projects/bio-bwa)**  
```
bwa mem refcontigs.fa trimmed.ind1.fq.gz > ind1.sam
```

3. **Convert .sam to .bam using [samtools](http://www.htslib.org)**  
```
samtools view -bS ind1.sam > ind1.unsorted.bam
```

4. **Sort .bam using samtools**    
```
samtools sort ind1.unsorted.bam > ind1.bam
```

5. **Add read group tags using [picard](https://github.com/broadinstitute/picard)**  
```
java -jar picard.jar AddOrReplaceReadGroups I=ind1.bam O=ind1.RG.bam SO=coordinate RGID=SeqRUN# RGLB=ind1.bam RGPL=illumina RGPU=ind1.bam RGSM=ind1.bam VALIDATION_STRINGENCY=LENIENT
```

6. **Index .bam file using samtools**  
```
samtools index ind1.RG.bam ind1.RG.bai
```

7. **Call variants using [GATK](https://software.broadinstitute.org/gatk/)**    

```
java -jar GenomeAnalysisTK.jar -R refcontigs.fa -T UnifiedGenotyper -I ind1.RG.bam -I ind2.RG.bam -I ind3.RG.bam -ploidy 2 -rf BadCigar -o data.vcf
```

# Filtering your VCF file

There are several existing tools for filtering VCF files (e.g., vcftools).  

Another option is to use **filterVCF.py**

This script filters a VCF file based on mapping quality (default minimum mapping quality is 20), missing data, and proportion of reads that match the Ref allele (across all individuals).

### Usage

> python filterVCF.py [input VCF] [output filtered vcf] [n indivs] [min indivs per site]

* input VCF: a VCF file output by GATK
* output filtered VCF: desired output VCF filename
* n indivs: number of individuals included in the VCF file
* min indivs per site: minimum number of individuals with data for site to be written to output VCF

### Example

```
python filterVCF.py toy.vcf toy.filtered.vcf 291 100
```


# Simulations

**data.simulations.py** simulates phenotype and genotype data for a locus that affects phenotype and shallow sequence data. 

### Usage

> python data.simulations.py [outfile]

Critical variables in the script to be modified:
* `nreps`: number of replicate simulations
* `nsamples`: number of individuals sampled from the population
* `meandepth`: mean sequence read depth
* `true_q`: population allele frequency of the simulated locus
* `true_a`: additive effect of the simulated locus (effect of Alt allele)
* `true_muRR`: mean phenotype of homozygous (Ref allele) individuals
* `Vp`: phenotypic variance in the population

All these variables will have to be changed, depending on your experimental design, details of phenotype, etc.  
Note that if you model a locus with no effect on phenotype, set `true_a` = 0 and `true_muRR` = population mean phenotype.  
After modeling phenotype and genotype across all individuals, as well as missing data across individuals, this script then uses the same likelihood model as in **MLE.phenoassoc.py** to test for an association between genotype and phenotype.


### Example

```
python data.simulations.py sim.out.txt
```

This file outputs the following data for each replicate simulation (MLE = maximum likelihood estimate):
* rep: index for replicate
* q: replicate's MLE for allele frequency for reference allele (*q*)
* muRR: replicate's MLE for mean phenotype for ref allele homozygote under hypothesis that genotypes differ in mean phenotype
* sigma: replicate's MLE for phenotypic variance
* a: replicate's MLE for additive effect of Alt allele
* LRT: likelihood ratio test statistic value


# Estimating linkage disequilibrium (LD) between sites on the same sequencing read

Scripts to estimate LD using single end Illumina sequence data using the LDx approach of [Feder et al. (2012)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048588) are called **pyLDx.phase1.py** and **pyLDx.phase2.py**

## Phase 1

In the first phase, we find all instances of sites contained in a given VCF that occur on the same sequencing read. This information is found in .sam files generated by a read alignment program such as BWA. *This phase generates a large file.*

### Usage

> python pyLDx.phase1.py [input vcf] [input sam] [output file]

## Phase 2

In the second phase, we compile haplotype information from the output of phase 1 to calculate LD between pairs of sites sequenced on the same SE reads.

### Usage

> python pyLDx.phase2.py [phase1 file] [output file]

### Example

```
python pyLDx.phase1.py input.vcf input.sam haplos.txt
python pyLDx.phase2.py haplos.txt LDresults.txt
```

This file outputs the following data for pair of sites in the vcf that were sequenced on the same read:
* contig: identifier of contig or scaffold in ref genome
* site1: position of first SNP in the pair
* site2: position of second SNP in the pair
* distance: distance in basepairs between site1 and site2
* freq_site1: proportion of cases where site1 shows Ref base -> (xRR + xRA)/(xRR + xRA + xAR + xAA)
* freq_site2: proportion of cases where site2 shows Ref base -> (xRR + xAR)/(xRR + xRA + xAR + xAA)
* xRR: count of cases where both site1 and site2 show the Ref base
* xRA: count of cases where site1 shows Ref base and site2 shows Alt base
* xAR: count of cases where site1 shows Alt base and site2 shows Ref base
* xAA: count of cases where both site1 and site2 show the Alt base
* r2: measure of LD, following the calculation in Feder et al. (2012)
