# CIBERSORT

Here I am checking the proportion of genes in the `LM22` signature matrix are included from the Twist Exome panel.
LM22 signature matrix was downloaded from the [CIBERSORT website](https://cibersortx.stanford.edu/download.php).

It contains the following 22 cell types:

* B cells naive
* B cells memory  
* Plasma cells
* T cells CD8
* T cells CD4 naive
* T cells CD4 memory resting
* T cells CD4 memory activated
* T cells follicular helper
* T cells regulatory (Tregs)
* T cells gamma delta
* NK cells resting
* NK cells activated
* Monocytes
* Macrophages M0
* Macrophages M1
* Macrophages M2
* Dendritic cells resting
* Dendritic cells activated
* Mast cells resting
* Mast cells activated 
* Eosinophils
* Neutrophls

The proportion of genes from the LM22 signature matrix found in the Twist Illumina Exome bed file are as follows: 

```console
Number of genes in LM22 signature found in Illumina Twist Exome kit: 510
Number of genes in LM22 signature not found in Illumina Twist Exome kit: 37
Percentage of genes in LM22 signature found in Illumina Twist Exome kit: 93.2
Genes in LM22 signature not found in file 2: {'HIST1H2AE', 'GPR97', 'HMGB3P30', 'FAIM3', 'SEPT8', 'KIRREL', 'LINC00597', 'CLCA3P', 'GUSBP11', 'HIST1H2BG', 'LOC100130100', 'EMR1', 'FLJ13197', 'FAM198B', 'LAT', 'LRMP', 'GSTT1', 'CXorf57', 'LOC126987', 'FAM212B', 'ZNF204P', 'KRT18P50', 'PCDHA5', 'KIAA0754', 'RPL3P7', 'SEPT5', 'LILRA3', 'LINC00921', 'FLT3LG', 'KIAA0226L', 'MARCH3', 'TARDBPP1', 'ATHL1', 'IGLL3P', 'FAM65B', 'EMR2', 'EMR3'}
```

According to the [CIBERSORT website tutorial](https://cibersortx.stanford.edu/tutorial.php):
```console
CIBERSORTx performs a feature selection and therefore typically does not use all genes in the signature matrix. It is generally ok if some genes are missing from the user’s mixture file. If less than 50% of signature matrix genes overlap, CIBERSORTx will issue a warning.
```

# Quantiseq

The signature matrix for [quantiseq](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0638-6) is the `TIL_10` signature matrix and has the following 10 cell types:

* B.cells
* Macrophages.M1
* Macrophages.M2
* Monocytes
* Neutrophils
* NK.cells
* T.cells.CD4
* T.cells.CD8
* Tregs
* Dendritic.cells

```console
Number of genes in TIL10 signature found in Illumina Twist Exome kit: 165
Number of genes in TIL10 signature not found in Illumina Twist Exome kit: 5
Genes in TIL10 signature not found in file 2: {'GUCY1A3', 'ADAM6', 'AATBC', 'AKAP2', 'FAM46C'}
Percentage of genes in TIL10 signature found in Illumina Twist Exome kit: 97.1
```

# xCell

[xCell](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1) contains a number of gene signatures on which gene signature scoring is carried out - it is not an absolute approach.

Here is the coverage:

```console
Number of genes in xcell signatures found in Illumina Twist Exome kit: 4936
Number of genes in xcell signatures not found in Illumina Twist Exome kit: 143
Percentage of genes in xcell signatures found in Illumina Twist Exome kit: 97.2
Genes in xcell signatures not found in file 2: {'KIAA1107', 'C2orf47', 'HOXA6', 'C16orf59', 'SEPT2', 'CECR5', 'SEPT7', 'PCDHA5', 'C14orf2', 'GPR27', 'ZNF154', 'TCL6', 'FAM134C', 'PLA2G16', 'KIAA1033', 'ZNF548', 'CNTD2', 'PCDHGA9', 'HIST1H2BL', 'HOXA4', 'KIAA0196', 'TSSC1', 'RARS', 'WHSC1', 'TSSK1B', 'PCDHGA11', 'C11orf57', 'LECT1', 'TMEM110', 'MKL1', 'RNF219', 'C14orf169', 'C21orf59', 'HIST1H3C', 'C6orf10', 'FAM21A', 'GARS', 'KIRREL', 'H3F3B', 'KIAA0101', 'NARFL', 'PCDHA10', 'PQLC2', 'C14orf1', 'IARS', 'ATP5G1', 'ATP5L', 'PNMAL1', 'FAM64A', 'H2AFY', 'EPRS', 'ATP5SL', 'ALOX12', 'MARS', 'H2AFZ', 'KARS', 'NUPL2', 'FAM65A', 'SEPT9', 'HIST1H2BB', 'C1orf61', 'MARCH8', 'C16orf58', 'C12orf49', 'C14orf105', 'MKL2', 'C19orf24', 'ATP5C1', 'SEPP1', 'CCDC53', 'ATP5G3', 'CTGF', 'TARS', 'WDR60', 'FAM175B', 'ATPIF1', 'PCDHGB5', 'FMO6P', 'C12orf10', 'FAM96B', 'CNPY2', 'CXorf36', 'COL4A3BP', 'NOV', 'KIAA0430', 'HIST1H3A', 'LIPT1', 'TCEB1', 'HARS', 'PARK2', 'C1orf123', 'ZCCHC11', 'C7orf43', 'HIST1H1D', 'WISP1', 'HIST1H2BM', 'T', 'C10orf76', 'PCDHA6', 'HIST1H1A', 'LOR', 'NARS', 'DFNA5', 'CXorf21', 'KIAA0368', 'IKBKAP', 'GLTSCR2', 'CPSF3L', 'HIST1H2BO', 'NAT6', 'UFD1L', 'MUM1', 'EGOT', 'AIM1', 'SHFM1', 'HOXA5', 'HIST1H4C', 'H2AFV', 'SDCCAG3', 'CD3EAP', 'HIST1H4F', 'C16orf45', 'ADSS', 'KIAA0125', 'C16orf62', 'H2AFX', 'C14orf166', 'FAM134B', 'TROVE2', 'ATP5J', 'HIST1H3D', 'C6orf25', 'C21orf2', 'SEPT10', 'CA7', 'ARSE', 'MIF', 'FAM105A', 'C1orf27', 'WARS', 'G6PC', 'MFSD7', 'DGCR14'}
```

# MCP Counter

MCP counter infers signature scores for the following cell types:

* T cells
* CD8 T cells
* Cytotoxic lymphocytes
* B lineage
* NK cells
* Monocytic lineage
* Myeloid dendritic cells
* Neutrophils
* Endothelial cells
* Fibroblasts

```console
Number of genes in mcp counter signature found in Illumina Twist Exome kit: 104
Number of genes in mcp counter signature not found in Illumina Twist Exome kit: 7
Percentage of genes in mcp counter signature found in Illumina Twist Exome kit: 93.7
Genes in mcp counter signature not found in file 2: {'KIR3DS1', 'MGC40069', 'CHRM3-AS2', 'FLT3LG', 'WFDC21P', 'C', 'KIR2DL3'}
```