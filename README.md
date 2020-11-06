This library provides a brief description of Gene-based test pipeline  
1. Quality Control(QC)  
2. Annotation  
3. Gene based test using SKAT function (SKAT R - package : <https://github.com/leeshawn/SKAT>)  

-------------------------------------------
# 1. Quality Control(QC) - vcftools  
   * Original filtered vcf file  
   
<pre>
<code>
--filter-name LowReadPosRankSum --filter-expression "ReadPosRankSum < -2.0"   
--filter-name LowMQRankSum --filter-expression "MQRankSum < -2.0"   
--filter-name LowQual --filter-expression "QUAL < 30.0"  
--filter-name QD --filter-expression "QD < 3.0"  
--filter-name FS --filter-expression "FS > 30.0"  
--filter-name MQ --filter-expression "MQ < 30.0"  
--filter-name DP --filter-expression "DP < 10"  
--genotype-filter-name DP --genotype-filter-expression "DP < 10"  
--genotype-filter-name GQ --genotype-filter-expression "GQ < 10.0"  
</code>
</pre>

   * Additional Quality Control  
        * Vcftools  
        Vcftools is a suite of functions for use on genetic variation data in the form of VCF and BCF files.  
        We can download the vcftools here : https://vcftools.github.io/downloads.html or    
        it can be installed using Anaconda.
        <pre>
        <code>
        conda activate vcf_env    # creating vcf virtual environment
        conda intall -c bioconda vcftools     # install vcftools
        </code>
        </pre>  


   * a) Hardy-Weinberg equilibrium(HWE)  
        * HWE < 10e-6  
        
   * b) Missing rate  
        * Missing rate < 0.1  

   
   


# 2. Annotation - Annovar
ANNOVAR (ANNOtate VARiation) is a bioinformatics software tool for the interpretation and prioritization of single nucleotide variants (SNVs), insertions, deletions, and copy number variants (CNVs) of a given genome. It has the ability to annotate human genomes hg18, hg19, and hg38.


![annovar](https://user-images.githubusercontent.com/73377376/97069199-74801580-1609-11eb-8775-0b07cadf878d.png)


### Example of basic workflow

1. Download annotation database

    Here we can see the available databases : <http://annovar.openbioinformatics.org/en/latest/user-guide/download/>  
    You can download the reference genome database that suits your study (Ex. hg19(GRCh37), hg38(GRCh38) .. )

   * Download the hg38 database
<pre>
<code>
        ./perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
        # Or if perl is already installed,
        ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
</code>
</pre>

2. Run annotate-varation.pl and proceed with the annotation work using the vcf file

<pre>
<code>
        ./perl table_annovar.pl  [VCF_FILENAME.vcf]  humandb/ --outfile [FILENAME] --buildver hg38 
        --protocol refGene --operation g –vcfinput
        # Or if perl is already installed,
        ./table_annovar.pl  [VCF_FILENAME.vcf]  humandb/ --outfile [FILENAME] --buildver hg38 
        --protocol refGene --operation g –vcfinput
</code>
</pre>  
* VCF_FILENAME.vcf : input vcf filename  
* FILENAME : name of the file you want to save  

For more information, please refer to this tutorial : <http://annovar.openbioinformatics.org/en/latest/user-guide/startup/#annotate_variationpl>

If the process above is successfully performed, the following file will be created. > [FILENAME].hg38_multianno.txt.  
The following analysis can be carried out using the multianno.txt file above.


-------------------------------------


# 3. Gene based test - SKAT

### SKAT_gene_SSD_All_functions.R

   The libraries used : SKAT, data.table, tibble, dplyr  
   **Before executing this function, place the plink or plink2 file on the current directory**  
   
   
    
  * Convert_Annovar function  
      * Converting the chromosome name from chr1-22, chrX-Y to 1-22, X-Y and read the annotation result as a csv  
        (the chromosome name may be inconsistent at times, so a conversion may be necessary)  
        
        *Input = (Annovar result multianno.txt file, Name to save after converting)*  
        
        <pre>
        <code>
        Convert_Annovar('SNUH_annotation.hg38_multianno.txt', 'SNUH_converted_multianno.csv')
        </code>
        </pre>  
        
  * Var_Info function
      * Using the annovar result file(the converted file above), create variant information table for each gene.  
        It shows the number of variants with different functions (e.g. exonic, splicing). Exonic function can be seen for the exonic variants.  
        Also, generate loss of function column by adding the results of five variants(frameshift deletion, frameshift insertion, 
        startloss, stopgain, stoploss)
        
        
        *Input = (Annovar result converted.txt, Name to save after converting)*  

        
        <pre>
        <code>
        Var_Info('SNUH_converted_multianno.csv', 'SNUH_var_info.txt')
        </code>
        </pre>  

        
  * SKAT_gene_SSD function
      * Gene-wise SKAT (Security Kernel Association Test) analysis to confirm the significance of genes  
        Before executing this function, you must convert your vcf file used in the annoar software into a bfile (bed, bim, fam) via plink.  
        For example, "plink --vcf [VCF_FILENAME].vcf --make-bed --out [BFILE_NAME]"
      
        *Input = (Annovar result converted.txt, bfile(bed, bim, fam) name, Name to save after processing, cov=NULL, method='SKATO', weights.beta=c(1,25), weights=NULL,  
                  Is.binary, genefunc=c(), exonicfunc=c(), number=1000, leaveSSD=F, plinkver=2)*    
        
        * **cov** : Name of covariate cov.file, default=NULL  
        * **method** : {SKAT, SKATO, Burden}. default='SKATO'   
        * **weights.beta** : The degree of weighting according to Minor Allele Frequency(MAF)  
        * **weights** : Weights in SKAT function.  
        * **Is.binary** : Check if phenotype is binary or not  
        * **genefunc** : Gene variables to be used for the SKAT test   
            * gene_varlist = 
                { downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_exonic;splicing, ncRNA_intronic, ncRNA_splicing, ncRNA_UTR5, splicing, upstream,                   upstream;downstream, UTR3, UTR5, UTR5;UTR3 }   
        * **exonicfunc** : Exonic variables to be used for the SKAT test   
            * exonic_varlist = 
                { frameshift deletion, frameshift insertion, nonframeshift deletion, nonframeshift insertion, nonsynonymous SNV, startloss, stopgain, stoploss, synonymous SNV } 
        * **number** : The number of genes to be used for each time. If n=1000, analyze 1000 genes at a time(Setting for computation). Default = 1000  
        * **leaveSSD** : If leaveSSD=F, save the results separately for each number. default=F  
                           To use the Oneset_genotype_SSD function below, you need to set leaveSSD=T
        * **plinkver** : plink version. default=2 (plink2)  
        
        
        <pre>
        <code>
        SKAT_gene_SSD('SNUH_converted_multianno.csv','SNUH','SNUH_result.txt', cov=NULL, method='SKATO', weights.beta=c(1,25),  
                       weights=NULL, Is.binary=T, genefunc=c('exonic','splicing','exonic.splicing'),exonicfunc=c('nonsynonymous SNV'), n=1000, leaveSSD=T)
        </code>
        </pre>  
                  
  * Oneset_genotype_SSD function
      * Print the genotype matrix of the relevant gene.  
        This function can be used only when the leaveSSD=T of SKAT_gene_SSD function.
          
        *Input = (gene name, Annovar result converted.txt, bfile(bed, bim, fam) name, number used in SKAT function)*
        
        <pre>
        <code>
        Oneset_Genotype_SSD('FYB2','SNUH_converted_multianno.csv','SNUH',1000)
        </code>
        </pre>  
 
  * SKAT_gene_SSD_specific function
      * This function is a modified version of the existing SKAT_gene_SSD function.  
        With the gene=c() parameter, a specific gene list can be inserted into the input to conduct group analysis only for these genes.  
    
        *Input = (Annovar result converted.txt, bfile(bed, bim, fam) name, **gene=c()**, Name to save after processing, cov=NULL, method='SKAT', weights.beta=c(1,25),  
                  weights=NULL, Is.binary, genefunc=c(), exonicfunc=c(), number=1000, leaveSSD=F, plinkver=2)*        
        
        * **gene** : Gene list to be used for gene group analysis 
        * **weights.beta** : The degree of weighting according to Minor Allele Frequency(MAF)
        * **weights** : Weights in SKAT function.  
        
        
        <pre>
        <code>
        SKAT_gene_SSD_specific('SNUH_converted_multianno.csv','SNUH', gene=c('HLA-DRB1', 'HLA-DRB5', 'HLA-B'), 'SNUH_result.txt', cov=NULL, method='SKATO',  
                                weights.beta=c(1,25), weights=NULL, Is.binary=T, genefunc=c('exonic','splicing','exonic.splicing'),exonicfunc=c('nonsynonymous SNV'), n=1000)
        </code>
        </pre>  
              
              
# Code examples of SNUH
   * Files and parameters  
        **environment** = Execution of R codes under linux environment  
        **multianno file** = SNUH_annotation.hg38_multianno.txt  
        **bfile name** = SNUH (SNUH.bed, SNUH.bim, SNUH.fam)  
        **covar** = NULL  
        **plink version** = plink2  
        **method** = SKATO  
        **gene functions** = c('exonic', 'exonic;splicing', 'ncRNA_exonic;splicing' , 'ncRNA_splicing' , 'splicing')  
        **exonic functions** = c('frameshift deletion', 'frameshift insertion', 'nonframeshift deletion', 'nonframeshift insertion',  
                           'nonsynonymous SNV', 'startloss', 'stopgain', 'stoploss')  
                             
                             
        <pre>
        <code>
     
        # Read_Annovar funciton  
        Convert_Annovar('SNUH_annotation.hg38_multianno.txt','SNUH_converted_multianno.csv') # Then, we will get the converted csv file  
        
        # SKAT_gene_SSD function  
        SKAT_gene_SSD('SNUH_converted_multianno.csv','SNUH','SNUH_example_result.txt',method='SKATO',weights.beta=c(1,25),Is.binary=T,  
                       genefunc=c('exonic', 'exonic;splicing', 'ncRNA_exonic;splicing' , 'ncRNA_splicing' , 'splicing'),  
                       exonicfunc=c('frameshift deletion', 'frameshift insertion', 'nonframeshift deletion', 'nonframeshift insertion',  
                       'nonsynonymous SNV', 'startloss', 'stopgain', 'stoploss'),leaveSSD=F)  
        # Then, we can check the results using the SNUH_example_result.txt file generated from the SKAT_gene_SSD function.  
        </code>
        </pre>  


          
## References
* Lee, S., Emond, M.J., ..., and Lin, X. (2012). Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. AJHG, 91, 224-237.  
* Wu, M., Lee, S., Cai, T., Li, Y., Boehnke, M. and Lin, X. (2011). Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). AJHG, 89, 82-93.  
* Lee, S., Fuchsberger, C., Kim, S., Scott, L. (2016) An efficient resampling method for calibrating single and gene-based rare variant association analysis in case-control studies, Biostatistics, 17, 1-15.
