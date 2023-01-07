---
title: "Phylogenetics with `phruta`"
author: "Cristian Román-Palacios"
date: "2023-01-06"
slug: "phylogenetics-with-phruta"
tags:
- biology
- biodiversity
- bioinformatics
- genbank
- phylogeny
package_version: 0.1.3
description: Using `phruta` to reconstruct the phylogeny of the new world Quail
---




<small><em>This post discusses the phylogenetic workflow implemented in the newly released `phruta` `R` package. 
As the project evolves, we will post updates to document features and technical details.
For more information, visit the [project website](https://ropensci.github.io/phruta/).</em></small>


## What is `phruta` and how can you use it?

The `phruta` package is primarily designed to simplify the basic phylogenetic computational workflow in `R`. `phruta` is expected to allow scientists to assemble molecular databases or phylogenies for particular taxonomic groups, with minimal complexity and maximal reproducibility. All the code in `phruta` runs within `R`. 


## Using `phruta` to assemble a phylogeny from scratch

Below, we will use `phruta` to assemble a phylogeny for the new world Quails. Species in this group are classified in the family Odontophoridae, a clade including nearly 34 extant species classified in 10 genera. Let's start by loading `phruta`:


```r
library(phruta)
```



So far, we have decided the taxonomic makeup of our analyses. We will now use `phruta` to figure out what genes are well sampled in GenBank for both the ingroup and outgroup:


```r
gs.seqs <- gene.sampling.retrieve(organism = c("Odontophoridae",
    "Ptilopachus", "Polyplectron"), speciesSampling = TRUE, npar = 6,
    nSearchesBatch = 500)
```

For the search terms, `phruta` was able to retrieve the names for 1 gene regions from GenBank. We will now generate a preliminary summary of the accession numbers retrieved for the combination of target taxa and gene regions. For simplicity, we will focus on analyzing gene regions that are sampled in \>20% of the species (`targetGenes` data.frame): 


```
##                           Gene Sampled in N species PercentOfSampledSpecies
## 1 NADH dehydrogenase subunit 2                   44                     100
```

Next, we will create a new object, `acc.table`, that will later be used to download the relevant gene sequences from GenBank. 


```r
targetGenes <- gs.seqs[gs.seqs$PercentOfSampledSpecies > 20,
    ]
acc.table <- acc.table.retrieve(clades = c("Odontophoridae",
    "Ptilopachus", "Polyplectron"), genes = targetGenes$Gene,
    speciesLevel = TRUE, npar = 6, nSearchesBatch = 500)
```

So far, we have only assembled a sampling matrix. Let's now download all the sequences listed in the accessions table generated above:


```r
sqs.downloaded <- sq.retrieve.indirect(acc.table = acc.table,
    download.sqs = FALSE)
```

Now, let's make sure that we are only including sequences that are reliable and from species that we are actually interested in analyzing. 



```r
sqs.curated <- sq.curate(filterTaxonomicCriteria = "[AZ]", kingdom = "animals",
    sqs.object = sqs.downloaded, removeOutliers = FALSE)
```

Running the `sq.curate()` function will create an object of class `list` (i.e. `sqs.curated`) that includes (1) the curated sequences with original names, (2) the curated sequences with species-level names (`renamed_*` prefix), (3) the accession numbers table (`AccessionTable`), and (4) a summary of taxonomic information for all the species sampled in the files. 


```
##                     OriginalNames     AccN                Species
## 1 MZ476322 Callipepla californica MZ476322 Callipepla_californica
## 2    MZ476314 Callipepla gambelii MZ476314    Callipepla_gambelii
## 3    EU166949 Colinus virginianus EU166949    Colinus_virginianus
## 4      AF222544 Colinus cristatus AF222544      Colinus_cristatus
## 5   KR732857 Colinus nigrogularis KR732857   Colinus_nigrogularis
## 6    KR732856 Dendrortyx barbatus KR732856    Dendrortyx_barbatus
##                           file
## 1 NADH dehydrogenase subunit 2
## 2 NADH dehydrogenase subunit 2
## 3 NADH dehydrogenase subunit 2
## 4 NADH dehydrogenase subunit 2
## 5 NADH dehydrogenase subunit 2
## 6 NADH dehydrogenase subunit 2
```

From here, we will align the sequences that we just curated using `sq.aln()` with default parameters:


```r
sqs.aln <- sq.aln(sqs.object = sqs.curated)
```

The masked alignments are shown below.



{{< figure src = "UF.Cur.png" alt = "Sequence alignment" class = "center" >}}

Phylogenetic inference in `phruta` is conducted using the `tree.raxml()` function. To use this function, we will have to export our sequence alignments locally. We will follow the same folder structure as if we were exporting everything locally. Specifically, our sequence alignments will located in `2.Alignments` and we will exclusively export the alignments that were masked:




We are now ready to run RAxML from `phruta`:


```r
outgroup <- sqs.curated$Taxonomy[sqs.curated$Taxonomy$genus ==
    "Polyplectron", ]

tree.raxml(folder = "2.Alignments", FilePatterns = "Masked_",
    raxml_exec = "raxmlHPC", Bootstrap = 100, outgroup = paste(outgroup$species_names,
        collapse = ","))
```

```
## [1] "raxmlHPC -T 4 -f a -p 1234 -x 1234 -m GTRGAMMA -o Polyplectron_inopinatum,Polyplectron_napoleonis,Polyplectron_chalcurum,Polyplectron_bicalcaratum,Polyplectron_malacense,Polyplectron_germaini,Polyplectron_katsumatae -k   -N 100   -s phruta.phy -n phruta"
```

The resulting phylogeny from these analyses is presented below.



{{< figure src = "raxml_ingroup.png" alt = "RAxML tree" class = "center" >}}


Finally, let's perform tree dating in our phylogeny using secondary calibrations extracted from [Scholl and Wiens (2016)](https://royalsocietypublishing.org/doi/pdf/10.1098/rspb.2016.1334). I am only using this study because it has a large phylogeny but I expect to replace it in the near future:


```r
dir.create("1.CuratedSequences")
write.csv(sqs.curated$Taxonomy, "1.CuratedSequences/1.Taxonomy.csv")

tree.dating(taxonomyFolder = "1.CuratedSequences", phylogenyFolder = "3.Phylogeny",
    scale = "treePL")
```

The resulting time-calibrated tree is presented below.




{{< figure src = "phylo_ingroup.png" alt = "Timetree" class = "center" >}}

This is the end of this brief tutorial! Hope that you've found this tutorial to be informative and useful. Please [get in touch](https://cromanpa94.github.io/cromanpa/) if you find any issues, would like to help me improve this tool, would like to see any additional feature, etc. 

## A `shiny` app for `phruta`!!

You can also use the same functions implemented in `phruta` in a graphical environment. Please visit the repo for `salphycon`, the associated `shiny` app, and enjoy a simplified version of `phruta`. You can find the repo [here](https://github.com/cromanpa94/salphycon).


## Acknowledgements

Anna Krystalli(https://github.com/annakrystalli), [Rayna Harris](https://github.com/raynamharris), and [Frederick Boehm](https://github.com/fboehm) provided excellent comments during the peer review process of `phruta` in [ROpenSci](https://github.com/ropensci/software-review/issues/458). The author thanks [Heidi E. Steiner](https://heidiesteiner.netlify.app/) for proofreading the vignettes and documentation in `phruta` in addition to early versions of this manuscript.






