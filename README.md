## ADAGEpath

ADAGEpath provides necessary functions to perform ADAGE-based signature analysis.

You can install the latest development version from github with
``` r
install.packages("devtools")
devtools::install_github("ADAGEpath")
```
### ADAGE introduction
ADAGE (or eADAGE) models are neural networks trained unsupervisedly on large
publically available gene expression compendia. The gene expression signatures
extracted by ADAGE models resemble human-annotated biological pathways and cover
a significant amount of existing pathways [1,2]. In addition to signatures well
characterized by known pathways, ADAGE also extracts signatures that may
bring novel biological insights.

Please refer to the following papers if you want to learn more about ADAGE:  
1. [ADAGE-Based Integration of Publicly Available Pseudomonas aeruginosa Gene
Expression Data with Denoising Autoencoders Illuminates Microbe-Host
Interactions](http://msystems.asm.org/content/1/1/e00025-15)  
2. [Unsupervised extraction of functional gene expression signatures in the
bacterial pathogen Pseudomonas aeruginosa with eADAGE](
http://biorxiv.org/content/early/2016/12/02/078659)

### Preloaded data
ADAGEpath currently only supports one organism *Pseudomonas aeruginosa* and is
preloaded with a *P.a.* expression compendium (`PAcompendium` and its probe
quantile distribution `probedistribution`) and an eADAGE
model (`eADAGEmodel`) trained on the compendium. The package is also preloaded
with *P.a.* gene (`geneinfo` and `PAO1orthologs`) and operon (`operons`)
information.

### ADAGE signature analysis
An ADAGE signature analysis usually includes the following steps:

**Data loading**  
ADAGEpath currently supports raw microarray data in CEL format and processed
microarray or RNAseq data. Use function `load_dataset()` to load your own
dataset or datasets publicly available in ArrayExpress.

Since ADAGE only accepts expression values in the (0,1) range, we linearly
transform expression values to be between 0 and 1 using the function
`zeroone_norm()`.

**Signature activity calculation**  
We next calculate each signature's activity for each sample in the dataset
with the function `calculate_activity()`.

**Active signature detection**  
Now we want to identify signatures whose activities strongly vary with sample
phenotypes, such as signatures whose activities being significantly different in
a comparison of interest. We recommend using the more robust
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) approach
to test differential activation when sample size is small.
To facilitate the most frequent two-group comparison,
we wrapped a two-group limma test into the function `build_limma()`. You can
visualize the limma two-group test results using `plot_volcano()` and get
active signatures from the limma test result using `get_active_signatures()`.
We also provide examples of analyzing time-course experiments or factorial-design
experiments using limma in the vignettes.  
You can also use other statistical tests to identify active signatures.  
`plot_activity_heatmap()` generates a heatmap showing how signature activity
changes across samples.

**Signature overlap examination**  
To be robust to noise, ADAGE models would sometimes construct signatures that
have overlapping genes. We can check whether the active signatures identified above
overlap with others using `plot_signature_overlap()`. If there is a group of
signatures largely overlap, to reduce the number of signatures to look at,
we can calculate the marginal activity of each signature using
`calculate_marginal_activity()`, which is the remaining
activity of a signature after removing overlapping genes of another signature.
We can visualize whether a signature is still strongly active after the impact
of another signature has been removed using `plot_marginal_activation()`. Examples
of signature overlap examination are in vignettes **ArrayExpress-example** and
**Time-course-example**.

**Signature interpretation and visualization**  
Finally, to get a detailed view of a signature or a group of signatures,
we can check what genes are contained in them using
`annotate_genes_in_signatures()` and visualize these genesf through a
gene-gene network using `visualize_gene_network()`. We can also download
existing KEGG pathways using `fetch_geneset()` and associate
signatures to known KEGG pathways using `annotate_signatures_with_genesets()`.

### Vignettes
The package comes with 5 vignettes introducing how ADAGE signature analysis
can be performed for different types of datasets.

- **ArrayExpress-example**: provides an example of loading in a dataset available
on ArrayExpress and performing a complete ADAGE signature analysis workflow on
a experimental design with two phenotype groups.

- **User-input-example**: provides an example of loading in a local dataset
and performing a complete ADAGE signature analysis workflow on
a experimental design with two phenotype groups.

- **RNAseq-example**: provides an example of loading in a processed RNAseq dataset
available on ArrayExpress.

- **Factorial-design-example**: provides an example of detecting differentially
active signatures with limma for an experiment with factorial design.

- **Time-course-example**: provides an example of analyzing experiments with
complex experimental design. The first part looks at signatures that show
largest activity ranges across samples. The second part builds a limma model to
detect signatures with differential temporal patterns between treatment
and control.

