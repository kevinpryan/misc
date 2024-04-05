Here I am checking the proportion of genes in the LM22 signature matrix are included from the Twist Exome panel.
LM22 signature matrix was downloaded from the [CIBERSORT website](https://cibersortx.stanford.edu/download.php).

Here are the outputs of `count_overlap_lm22_twist_panel.py`:

```
517
number of genes in lm22: 547
number of genes not in lm22: 30
percentage of lm22 genes in twist exome kit: 94.51553930530164
genes missing: ['ATHL1', 'CLCA3P', 'CXorf57', 'EMR1', 'EMR2', 'EMR3', 'FAIM3', 'FAM212B', 'FAM65B', 'FLJ13197', 'GPR97', 'GSTT1', 'HIST1H2AE', 'HIST1H2BG', 'HMGB3P30', 'IGLL3P', 'KIAA0226L', 'KIAA0754', 'KRT18P50', 'LILRA3', 'LINC00597', 'LINC00921', 'LOC100130100', 'LOC126987', 'LRMP', 'MARCH3', 'RPL3P7', 'SEPT8', 'TARDBPP1', 'ZNF204P']
```

According to the [CIBERSORT website tutorial](https://cibersortx.stanford.edu/tutorial.php):
```
CIBERSORTx performs a feature selection and therefore typically does not use all genes in the signature matrix. It is generally ok if some genes are missing from the user’s mixture file. If less than 50% of signature matrix genes overlap, CIBERSORTx will issue a warning.
```
