## ADAGEpath

ADAGEpath provides necessary functions to perform ADAGE-based signature analysis.

### ADAGE introduction
ADAGE (or eADAGE) models are neural networks trained in an unsupervised manner on
large publicly available gene expression compendia. ADAGE aims at building
essential features that can reconstruct a compendium. We derived gene expression
signatures from ADAGE neural network nodes and found that these signatures
resemble human-annotated biological pathways and cover many existing pathways
[1,2]. In addition to signatures that match known pathways, ADAGE also extracts
signatures that may group genes in novel data-driven ways.

Please refer to the following papers if you want to learn more about ADAGE:
1. [ADAGE-Based Integration of Publicly Available Pseudomonas aeruginosa Gene
Expression Data with Denoising Autoencoders Illuminates Microbe-Host
Interactions](https://doi.org/10.1128/mSystems.00025-15)
2. [Unsupervised extraction of functional gene expression signatures in the
bacterial pathogen Pseudomonas aeruginosa with eADAGE](
https://doi.org/10.1101/078659)

### Preloaded data
ADAGEpath currently supports one organism *Pseudomonas aeruginosa*. The package
is preloaded with a *P.a.* expression compendium (`PAcompendium` and its probe
quantile distribution `probedistribution`) and an eADAGE
model (`eADAGEmodel`) trained on the compendium. The package is also preloaded
with *P.a.* gene (`geneinfo` and `PAO1orthologs`) and operon (`operons`)
information.

### ADAGE signature analysis
Signatures are gene sets derived from an ADAGE model. They are formed because
their genes are expressed coordinately in some samples in the compendium.
An ADAGE signature analysis aims to identify signatures in which the expression
of constituent genes are altered by an experimental treatment. Such signatures
may represent biological processes that are perturbed by the treatment. ADAGE
signature analysis usually includes the following steps:

#### Data loading
ADAGEpath currently supports raw microarray data in CEL format and processed
microarray or RNAseq data. Use function `load_dataset()` to load your own
dataset or datasets publicly available in ArrayExpress.

Since ADAGE only accepts expression values in the (0,1) range, we linearly
transform expression values to be between 0 and 1 using the function
`zeroone_norm()`.

#### Signature activity calculation
We next calculate each signature's activity for each sample in the dataset
with the function `calculate_activity()`.

#### Active signature detection
We next identify signatures whose activities strongly vary with treatments,
such as signatures whose activities are significantly different in conditions
of interest. We recommend using
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
to test differential activation, particularly when sample size is small.
To facilitate the most frequent two-group comparison,
we wrapped a two-group limma test into the function `build_limma()`. You can
visualize the limma two-group test results using `plot_volcano()` and get
active signatures from the limma test result using `get_active_signatures()`.
We also provide examples of analyzing time-course experiments or factorial-design
experiments using limma in the vignettes.
You can also use other statistical tests to identify active signatures.
`plot_activity_heatmap()` generates a heatmap showing how signature activity
changes across samples.

#### Signature overlap examination
To be robust to noise, ADAGE models would sometimes construct signatures that
have overlapping genes. We can check whether the active signatures identified above
overlap with others using `plot_signature_overlap()`. If there is a group of
signatures that largely overlap, to reduce the number of signatures to look at,
we can calculate the marginal activity of each signature using
`calculate_marginal_activity()`, which is the remaining
activity of a signature after removing overlapping genes of another signature.
We can visualize whether a signature is still strongly active after the impact
of another signature has been removed using `plot_marginal_activation()`. Examples
of signature overlap examination are in vignettes **ArrayExpress-example** and
**Time-course-example**.

#### Signature interpretation and visualization
Finally, to get a detailed view of a signature or a group of signatures,
we can retrieve their constituent genes using `annotate_genes_in_signatures()`
and visualize these genes through a
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
active signatures with limma for an experiment with factorial design. Such
design typically has two treatment factors and measurements at four treatment
combinations.

- **Time-course-example**: provides an example of analyzing experiments with
complex experimental design. The first part looks at signatures that show the
largest activity ranges across samples. The second part builds a limma model to
detect signatures with differential temporal patterns between treatment
and control.


### Installation

You can install the latest development version from github with
``` r
install.packages("devtools")
devtools::install_github("greenelab/ADAGEpath")
```
If you want to build the vignettes, run
``` r
devtools::install_github("greenelab/ADAGEpath", build_vignettes = TRUE)
```
You can also install it via BioInstaller
``` r
library(BiocInstaller)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("greenelab/ADAGEpath")
```
For earlier verisons of R ( < 3.5), instead use
``` r
library(BiocInstaller)
biocLite("greenelab/ADAGEpath")
```

#### Potential problems during installation:

- make: gfortran-4.8: No such file or directory.

    `ADAGEpath` depends on the Bioconductor package
`impute` that requires gfortrain-4.8. Install `gfortrain-4.8` to solve this problem.
You can follow instructions here to install `gfortran-4.8` on Mac
https://stackoverflow.com/questions/23916219/os-x-package-installation-depends-on-gfortran-4-8

- Error validating server certificate for 'https://hedgehog.fhcrc.org:443'

    Please follow this thread to accept the certificate
https://github.com/hadley/devtools/issues/1401

- `BiocInstaller` package not available

    You can install it via
    ```r
    devtools::install_bioc("BiocInstaller")
    ```
    or
    ```r
    source("https://bioconductor.org/biocLite.R")
    biocLite("BiocInstaller")
    ```
