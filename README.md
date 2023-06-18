# BCH339N: Predicting Tumor Phenotypes from Gene Expression

<br>
<br>

![slide_01](slides/slides-01.png)
The objective of our research is to predict various tumor phenotypes from gene expression data.


![slide_02](slides/slides-02.png)
Our presentation today will start out by outlining our two project
objectives.

After that, we’ll visualize how our patient samples cluster, then perform a differential expression analysis and gene set enrichment analysis.  

Next, we’ll use random forest predictions to assign tumor phenotypes.

And lastly we’ll discuss the broader impacts of our findings.


![slide_03](slides/slides-03.png)
Let’s start out by outlining our objectives and explaining where we sourced our data.


![slide_04](slides/slides-04.png)
We started out by searching The Cancer Genome Atlas for transcriptomic data of tumors that we predicted would be biologically diverse.

Our analyses are based around gene expression data for breast cancer samples, skin melanoma samples, low grade glioma samples (which is a type of tumor that occurs in non-neuronal nervous system cells, like the spinal cord), and lastly, samples from mesothelioma tumors (which is usually caused by asbestos exposure).


![slide_05](slides/slides-05.png)
Our first objective was to determine which underlying biological processes characterize each of the four tumor types.

Our second objective was to determine whether gene expression profiles could be used to predict the identity of a tumor.


![slide_06](slides/slides-06.png)
I’ll now go through the steps that we took to address the first objective.
The first step in working with high dimensional data like transcriptomic profiles is to cluster the information into interpretable dimensions. For this, I used the UMAP method for dimensionality reduction.


![slide_07](slides/slides-07.png)
UMAP is a non-linear algorithm that projects high dimensional data into two dimensions. So for example, the gene expression matrix with 60k genes and 1200 samples was able to be simplified into two dimensions using this algorithm.


![slide_08](slides/slides-08.png)
Why is dimensionality reduction useful?
It allows us to verify that the overarching differences between patient sample data is caused by biological groups like tumor type, as opposed to other underlying variables, perhaps like age or sex.

By verifying that clusters correspond to biological groups, the design of our differential gene expression experiment is more robust.


![slide_09](slides/slides-09.png)
Here is the two dimensional projection of expression counts.
We can see that for the most part, samples cluster based on their tumor identity.
It’s important to note that the distances between clusters do not have any meaning using UMAP, since the algorithm is non-linear.


![slide_10](slides/slides-10.png)
Now that we validated that tumor identity is a justifiable way to group our samples, we performed a differential expression analysis followed by a gene set enrichment analysis.


![slide_11](slides/slides-11.png)
Within this portion of the analysis, the first goal was to determine which genes were over-expressed in each tumor type.

Four DESeq experiments were designed so that log2Fold changes of a particular tumor type were compared to the remaining samples.

Once DESeq provided log2Fold changes and test statistics for each gene, a gene set enrichment analysis was performed to determine which biological pathways were over-expressed in each tumor type.

We used the Hallmark set of 50 well defined biological pathways for this analysis.

Statistics from DESeq were used to rank genes by their importance in each of the biological pathways.


![slide_12](slides/slides-12.png)
Here are the gene set enrichment results for breast cancer tumors.

The y axis defines the pathways expressed in the data, and the x axis describes the normalized enrichment scores. Negative enrichment scores indicate that the pathway was over-expressed in breast cancer samples.

Our results support that the most over-expressed genes belong to estrogen response pathways.


![slide_13](slides/slides-13.png)
Based on our data, we characterize breast cancer samples by their expression of estrogen response pathways.

This is a reasonable pathway since estrogen is responsible for female sex characteristics.


![slide_14](slides/slides-14.png)
Next, we look at the gene set enrichment results for skin cancer tumors.

Again, negative enrichment scores indicate that the pathway was over-expressed in skin cancer samples.

Our results support that the most over-expressed genes belong to MYC pathways.


![slide_15](slides/slides-15.png)
Based on our data, we characterize skin cancer samples by their expression of MYC pathways, which are oncogenic transcription factors.


![slide_16](slides/slides-16.png)
Next, we look at the gene set enrichment results for low grade glioma tumors.

Again, negative enrichment scores indicate that the pathway was over-expressed in glioma samples.

Our results support that the most over-expressed genes belong to hedgehog signaling pathways.


![slide_17](slides/slides-17.png)
Based on our data, we characterize low grade glioma samples by their expression of hedgehog signaling pathways, which play important roles in stem cell regulation.

This is a reasonable pathway since low grade gliomas occur in the spinal cord, which houses lots of stem cells.


![slide_18](slides/slides-18.png)
Next, we look at the gene set enrichment results for mesothelioma tumors.

Again, negative enrichment scores indicate that the pathway was over-expressed in mesothelioma samples.

Our results support that the most over-expressed genes belong to interferon response pathways.


![slide_19](slides/slides-19.png)
Based on our data, we characterize mesothelioma samples by their expression of interferon response pathways, which play important roles in cell immune response.

This is a reasonable pathway since asbestos exposure is a frequent cause of mesothelioma.


![slide_20](slides/slides-20.png)
![slide_21](slides/slides-21.png)
![slide_22](slides/slides-22.png)
![slide_23](slides/slides-23.png)
![slide_24](slides/slides-24.png)
![slide_25](slides/slides-25.png)
![slide_26](slides/slides-26.png)
![slide_27](slides/slides-27.png)
![slide_28](slides/slides-28.png)
![slide_29](slides/slides-29.png)
![slide_30](slides/slides-30.png)
![slide_31](slides/slides-31.png)
![slide_32](slides/slides-32.png)
![slide_33](slides/slides-33.png)
![slide_34](slides/slides-34.png)