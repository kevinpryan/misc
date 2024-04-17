This list of CTA/self antigen genes is this list used in [LENS: Landscape of Effective Neoantigens Software](https://academic.oup.com/bioinformatics/article/39/6/btad322/7162685). The list comes from the Cancer Testis Antigen Database (http://www.cta.lncc.br/) 

Here is the overlap between the Illumina Twist Exome kit and CTA/self antigen genes:

```{console}
Number of unique CTA/self antigen genes used in LENS: 214
Number of CTA/self antigen genes found in Illumina Twist Exome kit: 206
Number of CTA/self antigen genes not found in Illumina Twist Exome kit: 8
CTA/self antigen genes not found in file 2: {'CTAGE5', 'FAM46D', 'CCDC36', 'VENTXP1', 'DSCR8', 'SSX6', 'ZNF645', 'SSX9'}
Percentage of CTA/self antigen genes found in Illumina Twist Exome kit: 96.3
```

To reproduce these results, run:

```{bash}
bash get_files.sh
python count_cta_self_genes_illumina_panel_2.py
bash remove_files.sh
```