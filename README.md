This library provides a brief use of annovar software and a description of the function of gene-based analysis using SKAT.  
(SKAT R - package : <https://github.com/leeshawn/SKAT>)

-------------------------------------------

# ANNOVAR
ANNOVAR (ANNOtate VARiation) is a bioinformatics software tool for the interpretation and prioritization of single nucleotide variants (SNVs), insertions, deletions, and copy number variants (CNVs) of a given genome. It has the ability to annotate human genomes hg18, hg19, hg38.


![annovar](https://user-images.githubusercontent.com/73377376/97069199-74801580-1609-11eb-8775-0b07cadf878d.png)


### Example of basic workflow

1. Download annotation database

    Here we can see the available databases : <http://annovar.openbioinformatics.org/en/latest/user-guide/download/>

   * download hg38 database
<pre>
<code>
        ./perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
        # Or if perl is already installed,
        ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
</code>
</pre>

2. Run annotate-varation.pl and proceed with the annotation work using vcf file

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

If the above process is successfully performed, a file will be created [FILENAME].hg38_multianno.txt.  
The following analysis can be carried out using the above multianno.txt file.


-------------------------------------


# Gene based analysis using SKAT

SKAT_gene_SSD.R
  * Read_Annovar function  
      * Converting chromosome name from chr1-22, chrX-Y to 1-22, X-Y and read the annotation result as table  
        (sometimes chromosome name may be inconsistent, so conversion was performed)  
        
        *Input = (Annovar result multianno.txt file, Name to save after converting)*  
        
        <pre>
        <code>
        Read_Annovar('example.hg38_multianno.txt', 'converted_multianno.txt')
        </code>
        </pre>  
        
  * Var_Info function
      * Using the annovar result file(above converted file), create variant information table for each gene.  
        It shows the number of variants with different function (e.g. exonic, splicing). Exonic function can be seen for the exonic variants.  
        
        *Input = (Annovar result converted.txt, Name to save after converting, Order=T)*  
        
        * **Order=T** : Generate loss of function column by adding five variants(frameshift deletion, frameshift insertion, 
            startloss, stopgain, stoploss) and sort them in descending order
        
        <pre>
        <code>
        Var_Info('converted_multianno.txt', 'var_info.txt', Order=T)
        </code>
        </pre>  

        
  * SKAT_gene_SSD function
      * Gene-wise SKAT (Security Kernel Association Test) analysis to confirm the significance of genes
      
        *Input = (Annovar result converted.txt, bfile(bed, bim, fam) name, Name to save after processing, cov=NULL, method='SKAT', Is.binary, genefunc=c(), exonicfunc=c(), number=1000, leaveSSD=F, plinkver=2)*    
        
        * **cov** : Name of covariate cov.file, default=NULL  
        * **method** : {SKAT, SKAT-O, Burden}. default='SKAT'   
        * **Is.binary** : Check if phenotype is binary or not  
        * **genefunc** : Gene variables to be used for the SKAT test   
            * gene_varlist = 
                { downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_exonic;splicing, ncRNA_intronic, ncRNA_splicing, ncRNA_UTR5, splicing, upstream, upstream;downstream, UTR3, UTR5, UTR5;UTR3 }   
        * **exonicfunc** : Exonic variables to be used for the SKAT test   
            * exonic_varlist = 
                { frameshift deletion, frameshift insertion, nonframeshift deletion, nonframeshift insertion, nonsynonymous SNV, startloss, stopgain, stoploss, synonymous SNV }  
        * **number** : The number of genes to be used for each time. if n=1000, analyze 1000 genes at a time. default = 1000  
        * **leaveSSD** : If leaveSSD=F, save the results separately for each number. default=F  
        * **plinkver** : plink version. default=2 (plink2)  
        
        <pre>
        <code>
        SKAT_gene_SSD('converted_multianno.txt','SNU','SNU_result.txt', cov=SNU_covar.cov, method='SKAT',Is.binary=T,
                      genefunc=c('exonic','splicing','exonic.splicing'),exonicfunc=c('nonsynonymous SNV'), n=1000, leaveSSD=T)
        </code>
        </pre>  
                  
  * Oneset_genotype_SSD function
      * Enter the gene you want, and save the comparison results with the reference database
          
        *Input = (gene name, Annovar result converted.txt, bfile(bed, bim, fam) name, number used in SKAT function)*

        <pre>
        <code>
        Oneset_Genotype_SSD('FYB2','SNU_anno.txt','SNU',1000)
        </code>
        </pre>  
        
        
## Referenecs
* Lee, S., Emond, M.J., ..., and Lin, X. (2012). Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. AJHG, 91, 224-237.  
* Wu, M., Lee, S., Cai, T., Li, Y., Boehnke, M. and Lin, X. (2011). Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). AJHG, 89, 82-93.  
* Lee, S., Fuchsberger, C., Kim, S., Scott, L. (2016) An efficient resampling method for calibrating single and gene-based rare variant association analysis in case-control studies, Biostatistics, 17, 1-15.
