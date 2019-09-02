# Workflow Detail

* 'HES3\_GFP\_ESC': After raw data preprocessing

----------------------------------------------------------------
## 0. Raw Data Preprocess
#### raw data: ../data/raw/
* script: raw\_preprocess.R
1. Drop genes/enhancers with >75% zero across time points
2. Log2-transformation
3. Upper quartile normalization of sample
```bash
# usage: Rscript raw_preprocess.R [raw data filename] [target celltype name]
Rscript raw_preprocess.R HES3_GFP_ESC_raw_tpm HES3_GFP_ESC
```

----------------------------------------------------------------
## I. Data Preprocessing

* Command: 
```bash
python3.5 Preprocess.py [cell type]
```

### Input Data(Raw data): 
* CAGE FANTOM5 data (e.g. HES3\_GFP\_ESC)
### Preprocessing:
1. Select enhancers & genes.
2. Stitch enhancers (Currently Disabled).
3. Select Super-enhancers (by ROSE).
4. Map genes to super-enhancers.
5. Gene symbols are converted to HGNC symbols and only retain the first one.

### Data Produced: 
File Name | Detail
-------------------|-----------------------------
**Expression Profiles** |------------------------------------------------------
enhancers\_[cell type].tsv | Enhancer expression profile
genes\_[cell type].tsv | Gene expression profile
stitched\_enhancer.tsv | Stitched enhancer expression profile
super-enhancer.tsv | seRNA expression profile 
super-enhancer(\_avg)\_exp.tsv | Super-enhancers expression profile (avg means averaging replicates.)
**Mapping Tables** |------------------------------------------------------
enhancer\_loci.bed | For mapping genes to stitched Enhancers 
gene\_loci.bed | For stitching enhancers 
se2gene.bed | Super-enhancer to genes mapping table
symbol\_[cell type].tsv | All gene loci to gene symbol mapping table.
**Other** |------------------------------------------------------
stitched\_forcall.tsv | 'stitched\_enhancer.tsv' converted for ROSE. 

------------------------------------------------------------------
## II. Non-negative Matrix Factorization

* Command: 
```bash
python3.5 NMF.py [cell type] [average flag]
```

### Input Data: 
1. seRNA expression profile (super-enhancer.tsv)
2. genes expression profile (genes\_[celltype].tsv)

### NMF:
1. If avg\_flag then average replicate expression
2. State to sample:
    * Identify hidden states
    * (Disabled) determine stages.
   Use Agglomerative Clustering with average linkage from [scikit-learn](http://scikit-learn.org/0.18/modules/generated/sklearn.cluster.AgglomerativeClustering.html#sklearn.cluster.AgglomerativeClustering).        


### Data Produced:
File Name | Detail
----------|----------
(Disabled) stage[n]\_genes.tsv | Gene expression profile in stage 'n' (determined by NMF and Clustering).
state2seRNA.csv | States of each seRNA loci.
sample2state.csv | States of each sample.

------------------------------------------------------------------
## III. Co-expression Network (State-specific seRNA Selection)

* Command:
```bash
python3.5 fromW_coxnet.py [celltype] [# of states]
```

### Input Data: 
##### 1. Averaged gene expression profile
##### 2. Averaged super-enhancers expression profile
### fromW\_coxnet:
1. Select state-specific seRNA based on W matrix (state2seRNA.csv).
2. Filter out low-rank peak annotated genes.
3. Construct bipartite co-expression network between seRNAs & mRNAs.
4. Select state-specific mRNAs for each state-specific seRNA based on mutual rank. 

## Data Produced:
File Name | Detail
----------|----------
selected\_seRNA\_state[n].csv | state[n]-specific seRNA exp profile.
net\_corr\_mat[n].csv | Correlation matrix of seRNA & mRNA
state[n]-specific\_mRNA.csv | selected by mutual rank, denoted as highly co-expressed mRNA.
coexpress\_mrna\_state[n].pkl | dictionary: {seRNA: [mRNA gene symbols]}
coexpress\_locus\_state[n].pkl | dictionary: {seRNA: [mRNA loci]}

-------------------------------------------------------------------------
## IV-1. TF binding site enrichment of state-specific Super-Enhancer
* Note: Currently only for enrichment of 2 states
* Command:
```bash
python3.5 enrich_TF.py [celltype]
```
### Input Data:


1. wgEncodeRegTfbsClusteredV3.bed
2. super-enhancer_avg_exp.tsv
3. selected_seRNA_state{state_id}.csv


### Detail:


* Fisher's exact test to test TF binding site

### Data Produced:
File Name | Detail
----------|----------
TF_enrichment.png | Joint Plot of enrichment result 

----------------------------------------
## IV-2. Functional Enrichment Analysis
* Command:
```bash
python3.5 GO_SemSim.py --help
```

### Input Data: 


* Gene symbols of each module

### FEA:
1. Utilize topGO package to perform GSEA

-------------------------------------------------------
## IV-3. TF co-regulatory score enrichment of state-specific Super-Enhancer
* Command:
```bash
python3.5 coreg_TF.py [celltype]
```
### Input Data:


1. coexpress_locus_state{}.pkl
2. wgEncodeRegTfbsClusteredV3.bed
3. genes\_{}\_avg.tsv
### Detail:


1. Compute mean(Jaccard, Overlap_coeff) as co-regulatory score
2. Fisher's exact test 

### Data Produced:
File Name | Detail
----------|----------
coreg{}.pkl | dict: {seRNA: coreg_score, p-value}


----------------------------------------------