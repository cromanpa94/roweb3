---
title: Phylogenetics With `phruta`
author: Cristian Román Palacios
date: '2023-01-06'
slug: phylogenetics-with-phruta
tags:
  - Software Peer Review
  - packages
  - biology
  - biodiversity
  - bioinformatics
  - genbank
  - phylogeny
  - pipelines
  - workflow-automation
  - workflows
package_version: 0.1.3
description: Using `phruta` to reconstruct the phylogeny of the new world Quail
---
<script src="{{< blogdown/postref >}}index_files/kePrint/kePrint.js"></script>
<link href="{{< blogdown/postref >}}index_files/lightable/lightable.css" rel="stylesheet" />
<link href="{{< blogdown/postref >}}index_files/bsTable/bootstrapTable.min.css" rel="stylesheet" />
<script src="{{< blogdown/postref >}}index_files/bsTable/bootstrapTable.js"></script>
<script src="{{< blogdown/postref >}}index_files/kePrint/kePrint.js"></script>
<link href="{{< blogdown/postref >}}index_files/lightable/lightable.css" rel="stylesheet" />
<link href="{{< blogdown/postref >}}index_files/bsTable/bootstrapTable.min.css" rel="stylesheet" />
<script src="{{< blogdown/postref >}}index_files/bsTable/bootstrapTable.js"></script>




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
    "Ptilopachus", "Polyplectron"), npar = 6, nSearchesBatch = 500)
```

For the search terms, `phruta` was able to retrieve the names for 32 gene regions from GenBank. We will now generate a preliminary summary of the accession numbers retrieved for the combination of target taxa and gene regions. For simplicity, we will focus on analyzing the top 6 gene regions: 

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:500px; overflow-x: scroll; width:700px; "><table class="table table-striped table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Gene </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Sampled in N species </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PercentOfSampledSpecies </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cytochrome b </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 52.941176 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 47.058824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NADH-subunit 2 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 41.176471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> calbindin </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 35.294118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 29.411765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> high mobility group 17 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 29.411765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> transforming growth factor beta 2 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 17.647059 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> brain-derived neurotrophic factor </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 11.764706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cytochrome oxidase subunit I </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 11.764706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> eukaryotic translation elongation factor 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 11.764706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NADH dehydrogenase subunit 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 11.764706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neurotrophin 3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 11.764706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18S ribosomal RNA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> aldolase B fructose-bisphosphate </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> beta-nerve growth factor </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> calbindin 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> clathrin heavy polypeptide Hc </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crystallin alpha A </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cytochrome B </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cytochrome p450 family 19 subfamily A polypeptide 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> early growth response 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> eukaryotic elongation factor 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fibrinogen beta chain </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> muscle skeletal receptor tyrosine kinase </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> myoglobin </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nerve growth factor beta polypeptide </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pterin-4 alpha-carbinolamine dehydratase </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rhodopsin </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> serine/arginine-rich splicing factor 3 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tropomyosin 1 alpha </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> v-myc myelocytomatosis viral oncogene-like protein </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zinc finger protein </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.882353 </td>
  </tr>
</tbody>
</table></div>

Next, we will create a new object, `acc.table`, that will later be used to download the relevant gene sequences from GenBank. 

```r 
acc.table <- acc.table.retrieve(clades = c("Odontophoridae",
    "Ptilopachus", "Polyplectron"), genes = head(gs.seqs)$Gene,
    speciesLevel = TRUE, npar = 6, nSearchesBatch = 500)
```

So far, we have only assembled a sampling matrix. Let's now download all the sequences listed in the accessions table generated above:

```r 
sqs.downloaded <- sq.retrieve.indirect(acc.table = acc.table)
```

Now, let's make sure that we are only including sequences that are reliable and from species that we are actually interested in analyzing. 

```r 
sqs.curated <- sq.curate(filterTaxonomicCriteria = "[AZ]", kingdom = "animals",
    sqs.object = sqs.downloaded, removeOutliers = FALSE)
```

Running the `sq.curate()` function will create an object of class `list` (i.e. `sqs.curated`) that includes (1) the curated sequences with original names, (2) the curated sequences with species-level names (`renamed_*` prefix), (3) the accession numbers table (`AccessionTable`), and (4) a summary of taxonomic information for all the species sampled in the files. 

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:500px; overflow-x: scroll; width:700px; "><table class="table table-striped table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> OriginalNames </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> AccN </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Species </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> file </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> DQ485889 Callipepla gambelii </td>
   <td style="text-align:left;"> DQ485889 </td>
   <td style="text-align:left;"> Callipepla_gambelii </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AY952697 Colinus virginianus </td>
   <td style="text-align:left;"> AY952697 </td>
   <td style="text-align:left;"> Colinus_virginianus </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF068192 Cyrtonyx montezumae </td>
   <td style="text-align:left;"> AF068192 </td>
   <td style="text-align:left;"> Cyrtonyx_montezumae </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF252860 Oreortyx pictus </td>
   <td style="text-align:left;"> AF252860 </td>
   <td style="text-align:left;"> Oreortyx_pictus </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC495453 Callipepla squamata </td>
   <td style="text-align:left;"> KC495453 </td>
   <td style="text-align:left;"> Callipepla_squamata </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC857649 Colinus cristatus </td>
   <td style="text-align:left;"> KC857649 </td>
   <td style="text-align:left;"> Colinus_cristatus </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF028772 Callipepla californica </td>
   <td style="text-align:left;"> AF028772 </td>
   <td style="text-align:left;"> Callipepla_californica </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF028751 Callipepla douglasii </td>
   <td style="text-align:left;"> AF028751 </td>
   <td style="text-align:left;"> Callipepla_douglasii </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FR694138 Odontophorus speciosus </td>
   <td style="text-align:left;"> FR694138 </td>
   <td style="text-align:left;"> Odontophorus_speciosus </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FR694137 Odontophorus gujanensis </td>
   <td style="text-align:left;"> FR694137 </td>
   <td style="text-align:left;"> Odontophorus_gujanensis </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FR694136 Odontophorus capueira </td>
   <td style="text-align:left;"> FR694136 </td>
   <td style="text-align:left;"> Odontophorus_capueira </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AM236885 Ptilopachus nahani </td>
   <td style="text-align:left;"> AM236885 </td>
   <td style="text-align:left;"> Francolinus_nahani </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF330064 Polyplectron inopinatum </td>
   <td style="text-align:left;"> AF330064 </td>
   <td style="text-align:left;"> Polyplectron_inopinatum </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF330063 Polyplectron germaini </td>
   <td style="text-align:left;"> AF330063 </td>
   <td style="text-align:left;"> Polyplectron_germaini </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF330062 Polyplectron emphanum </td>
   <td style="text-align:left;"> AF330062 </td>
   <td style="text-align:left;"> Polyplectron_napoleonis </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AF330061 Polyplectron chalcurum </td>
   <td style="text-align:left;"> AF330061 </td>
   <td style="text-align:left;"> Polyplectron_chalcurum </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KX988022 Polyplectron bicalcaratum </td>
   <td style="text-align:left;"> KX988022 </td>
   <td style="text-align:left;"> Polyplectron_bicalcaratum </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EU839455 Polyplectron katsumatae </td>
   <td style="text-align:left;"> EU839455 </td>
   <td style="text-align:left;"> Polyplectron_katsumatae </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EU005262 Polyplectron schleiermacheri </td>
   <td style="text-align:left;"> EU005262 </td>
   <td style="text-align:left;"> Polyplectron_schleiermacheri </td>
   <td style="text-align:left;"> cytochrome b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732830 Cyrtonyx montezumae </td>
   <td style="text-align:left;"> KR732830 </td>
   <td style="text-align:left;"> Cyrtonyx_montezumae </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732829 Dactylortyx thoracicus </td>
   <td style="text-align:left;"> KR732829 </td>
   <td style="text-align:left;"> Dactylortyx_thoracicus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732828 Oreortyx pictus </td>
   <td style="text-align:left;"> KR732828 </td>
   <td style="text-align:left;"> Oreortyx_pictus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732827 Odontophorus erythrops </td>
   <td style="text-align:left;"> KR732827 </td>
   <td style="text-align:left;"> Odontophorus_erythrops </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732826 Odontophorus gujanensis </td>
   <td style="text-align:left;"> KR732826 </td>
   <td style="text-align:left;"> Odontophorus_gujanensis </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732825 Odontophorus stellatus </td>
   <td style="text-align:left;"> KR732825 </td>
   <td style="text-align:left;"> Odontophorus_stellatus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732824 Odontophorus capueira </td>
   <td style="text-align:left;"> KR732824 </td>
   <td style="text-align:left;"> Odontophorus_capueira </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732823 Odontophorus speciosus </td>
   <td style="text-align:left;"> KR732823 </td>
   <td style="text-align:left;"> Odontophorus_speciosus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732822 Odontophorus leucolaemus </td>
   <td style="text-align:left;"> KR732822 </td>
   <td style="text-align:left;"> Odontophorus_leucolaemus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732821 Odontophorus balliviani </td>
   <td style="text-align:left;"> KR732821 </td>
   <td style="text-align:left;"> Odontophorus_balliviani </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732820 Dendrortyx macroura </td>
   <td style="text-align:left;"> KR732820 </td>
   <td style="text-align:left;"> Dendrortyx_macroura </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732819 Philortyx fasciatus </td>
   <td style="text-align:left;"> KR732819 </td>
   <td style="text-align:left;"> Philortyx_fasciatus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732818 Callipepla gambelii </td>
   <td style="text-align:left;"> KR732818 </td>
   <td style="text-align:left;"> Callipepla_gambelii </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732817 Callipepla californica </td>
   <td style="text-align:left;"> KR732817 </td>
   <td style="text-align:left;"> Callipepla_californica </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732816 Callipepla douglasii </td>
   <td style="text-align:left;"> KR732816 </td>
   <td style="text-align:left;"> Callipepla_douglasii </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732815 Callipepla squamata </td>
   <td style="text-align:left;"> KR732815 </td>
   <td style="text-align:left;"> Callipepla_squamata </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732814 Colinus cristatus </td>
   <td style="text-align:left;"> KR732814 </td>
   <td style="text-align:left;"> Colinus_cristatus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732813 Colinus nigrogularis </td>
   <td style="text-align:left;"> KR732813 </td>
   <td style="text-align:left;"> Colinus_nigrogularis </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732812 Colinus virginianus </td>
   <td style="text-align:left;"> KR732812 </td>
   <td style="text-align:left;"> Colinus_virginianus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KR732832 Ptilopachus petrosus </td>
   <td style="text-align:left;"> KR732832 </td>
   <td style="text-align:left;"> Ptilopachus_petrosus </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749467 Polyplectron inopinatum </td>
   <td style="text-align:left;"> KC749467 </td>
   <td style="text-align:left;"> Polyplectron_inopinatum </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749466 Polyplectron germaini </td>
   <td style="text-align:left;"> KC749466 </td>
   <td style="text-align:left;"> Polyplectron_germaini </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749465 Polyplectron napoleonis </td>
   <td style="text-align:left;"> KC749465 </td>
   <td style="text-align:left;"> Polyplectron_napoleonis </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749464 Polyplectron chalcurum </td>
   <td style="text-align:left;"> KC749464 </td>
   <td style="text-align:left;"> Polyplectron_chalcurum </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749463 Polyplectron bicalcaratum </td>
   <td style="text-align:left;"> KC749463 </td>
   <td style="text-align:left;"> Polyplectron_bicalcaratum </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC778974 Polyplectron katsumatae </td>
   <td style="text-align:left;"> KC778974 </td>
   <td style="text-align:left;"> Polyplectron_katsumatae </td>
   <td style="text-align:left;"> 12S ribosomal RNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AY952687 Cyrtonyx montezumae </td>
   <td style="text-align:left;"> AY952687 </td>
   <td style="text-align:left;"> Cyrtonyx_montezumae </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AY952686 Colinus virginianus </td>
   <td style="text-align:left;"> AY952686 </td>
   <td style="text-align:left;"> Colinus_virginianus </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KF298498 Colinus cristatus </td>
   <td style="text-align:left;"> KF298498 </td>
   <td style="text-align:left;"> Colinus_cristatus </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749518 Polyplectron inopinatum </td>
   <td style="text-align:left;"> KC749518 </td>
   <td style="text-align:left;"> Polyplectron_inopinatum </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749517 Polyplectron germaini </td>
   <td style="text-align:left;"> KC749517 </td>
   <td style="text-align:left;"> Polyplectron_germaini </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749516 Polyplectron napoleonis </td>
   <td style="text-align:left;"> KC749516 </td>
   <td style="text-align:left;"> Polyplectron_napoleonis </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749515 Polyplectron chalcurum </td>
   <td style="text-align:left;"> KC749515 </td>
   <td style="text-align:left;"> Polyplectron_chalcurum </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749514 Polyplectron bicalcaratum </td>
   <td style="text-align:left;"> KC749514 </td>
   <td style="text-align:left;"> Polyplectron_bicalcaratum </td>
   <td style="text-align:left;"> calbindin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MT456615 Callipepla squamata </td>
   <td style="text-align:left;"> MT456615 </td>
   <td style="text-align:left;"> Callipepla_squamata </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MN986948 Callipepla californica </td>
   <td style="text-align:left;"> MN986948 </td>
   <td style="text-align:left;"> Callipepla_californica </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JQ175605 Odontophorus leucolaemus </td>
   <td style="text-align:left;"> JQ175605 </td>
   <td style="text-align:left;"> Odontophorus_leucolaemus </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JQ175602 Odontophorus gujanensis </td>
   <td style="text-align:left;"> JQ175602 </td>
   <td style="text-align:left;"> Odontophorus_gujanensis </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JQ174494 Colinus cristatus </td>
   <td style="text-align:left;"> JQ174494 </td>
   <td style="text-align:left;"> Colinus_cristatus </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AY666347 Colinus virginianus </td>
   <td style="text-align:left;"> AY666347 </td>
   <td style="text-align:left;"> Colinus_virginianus </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DQ433856 Oreortyx pictus </td>
   <td style="text-align:left;"> DQ433856 </td>
   <td style="text-align:left;"> Oreortyx_pictus </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DQ433416 Callipepla gambelii </td>
   <td style="text-align:left;"> DQ433416 </td>
   <td style="text-align:left;"> Callipepla_gambelii </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KM896482 Odontophorus capueira </td>
   <td style="text-align:left;"> KM896482 </td>
   <td style="text-align:left;"> Odontophorus_capueira </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JN801309 Cyrtonyx montezumae </td>
   <td style="text-align:left;"> JN801309 </td>
   <td style="text-align:left;"> Cyrtonyx_montezumae </td>
   <td style="text-align:left;"> cytochrome oxidase subunit 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AY952738 Cyrtonyx montezumae </td>
   <td style="text-align:left;"> AY952738 </td>
   <td style="text-align:left;"> Cyrtonyx_montezumae </td>
   <td style="text-align:left;"> high mobility group 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AY952737 Colinus virginianus </td>
   <td style="text-align:left;"> AY952737 </td>
   <td style="text-align:left;"> Colinus_virginianus </td>
   <td style="text-align:left;"> high mobility group 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749806 Polyplectron inopinatum </td>
   <td style="text-align:left;"> KC749806 </td>
   <td style="text-align:left;"> Polyplectron_inopinatum </td>
   <td style="text-align:left;"> high mobility group 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749805 Polyplectron germaini </td>
   <td style="text-align:left;"> KC749805 </td>
   <td style="text-align:left;"> Polyplectron_germaini </td>
   <td style="text-align:left;"> high mobility group 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749804 Polyplectron chalcurum </td>
   <td style="text-align:left;"> KC749804 </td>
   <td style="text-align:left;"> Polyplectron_chalcurum </td>
   <td style="text-align:left;"> high mobility group 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KC749803 Polyplectron bicalcaratum </td>
   <td style="text-align:left;"> KC749803 </td>
   <td style="text-align:left;"> Polyplectron_bicalcaratum </td>
   <td style="text-align:left;"> high mobility group 17 </td>
  </tr>
</tbody>
</table></div>

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
[1] "raxmlHPC -T 4 -f a -p 1234 -x 1234 -m GTRGAMMA -o Polyplectron_inopinatum,Polyplectron_germaini,Polyplectron_napoleonis,Polyplectron_chalcurum,Polyplectron_bicalcaratum,Polyplectron_katsumatae,Polyplectron_schleiermacheri -k   -N 100   -s phruta.phy -n phruta"
```

The resulting phylogeny from these analyses is presented below.



{{< figure src = "raxml_ingroup.png" alt = "RAxML tree" class = "center" >}}


Finally, let's perform tree dating in our phylogeny using secondary calibrations extracted from a large time-calibrated phylogeny with overlapping nodes:

```r 
dir.create("1.CuratedSequences")
write.csv(sqs.curated$Taxonomy, "1.CuratedSequences/1.Taxonomy.csv")

tree.dating(taxonomyFolder = "1.CuratedSequences", phylogenyFolder = "3.Phylogeny",
    scale = "treePL")
```

The resulting time-calibrated tree is presented below.




{{< figure src = "phylo_ingroup.png" alt = "Timetree for the selected taxa" class = "center" >}}


This is the end of this brief tutorial! Hope that you've found this tutorial to be informative and useful. Please [get in touch](https://cromanpa94.github.io/cromanpa/) if you find any issues, would like to help me improve this tool, would like to see any additional feature, etc. 

## A `shiny` app for `phruta`!!

You can also use the same functions implemented in `phruta` in a graphical environment. Please visit the repo for `salphycon`, the associated `shiny` app, and enjoy a simplified version of `phruta`. You can find the [repo here](https://github.com/cromanpa94/salphycon).


## Acknowledgements

[Anna Krystalli](https://github.com/annakrystalli), [Rayna Harris](https://github.com/raynamharris), and [Frederick Boehm](https://github.com/fboehm) provided excellent comments during the peer review process of `phruta` in [ROpenSci](https://github.com/ropensci/software-review/issues/458). The author thanks [Heidi E. Steiner](https://heidiesteiner.netlify.app/) for proofreading the documentation in `phruta`.






