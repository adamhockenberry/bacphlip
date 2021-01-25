---
title: 'BACPHLIP: Predicting bacteriophage lifestyle from conserved protein domains'
tags:
  - Python
  - bioinformatics
  - genomics
  - bacteriophage
  - microbiology
authors:
  - name: Adam J. Hockenberry^[Corresponding author]
    orcid: 0000-0001-9476-0104
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Claus O. Wilke
    affiliation: 1
affiliations:
 - name: Department of Integrative Biology, The University of Texas at Austin
   index: 1
date: 25 January 2021
bibliography: paper.bib

---

# Summary
Bacteriophages are broadly classified into two distinct lifestyles: temperate
phages are capable of a latent phase of infection within a host cell, whereas
virulent phages directly replicate and lyse host cells upon infection. Here, we
present `BACPHLIP`, a computational tool for predicting bacteriophage lifestyle
based solely on genome sequence data. On an independent test set of 423 phage
genomes, `BACPHLIP` achieves a classification accuracy of 98%, exceeding that of
the best previously available software (79%). 

**Availability and Implementation:** `BACPHLIP` is freely available
(https://github.com/adamhockenberry/bacphlip), with development code provided
separately (https://github.com/adamhockenberry/bacphlip-model-dev).

**Supporting Information:** Supporting text regarding biological findings and
full accuracy assessment is available in the model repository.

# Statement of need

Bacteriophages play important ecological roles
[@paez-espino_uncovering_2016; @nishimura_environmental_2017; @emerson_host-linked_2018; @daly_viruses_2019; @gregory_marine_2019],
influence both the physiology and evolution of host species
[@touchon_genetic_2016; @forterre_manipulation_2011; @carey_phage_2019],
and possess a number of unique traits [@ofir_contemporary_2018] that
are relevant to ongoing biotechnological and medical applications
[@dedrick_engineered_2019; @rodriguez-gonzalez_quantitative_2020]. The
availability of phage genome sequences has expanded in recent years due
to the development of meta-genomics
[@deng_viral_2014; @pope_whole_2015; @roux_viral_2015; @simmonds_virus_2017; @tisza_discovery_2020].

A particularly important phage phenotype is whether the phage has a
temperate (lysogenic) or virulent (lytic) lifestyle
[@bobay_adaptation_2013]. Phages are extraordinarily diverse and this
dichotomy is an over-simplification
[@abedon_bacteriophage_2008; @dion_phage_2020]; nevertheless lifestyle
classification is broadly recognized and important in numerous contexts
[@mavrich_bacteriophage_2017].

While phage lifestyle ultimately should be determined experimentally,
this is impractical for the entirety of newly discovered phage genomes.
Using computational methods, @mcnair_phacts_2012 developed a random
forest classifier (`PHACTS`) to predict lifestyle based off of sequence
similarity to a set of query proteins (randomly selected from the phage
proteome training set). More recently, @mavrich_bacteriophage_2017
described a computational classification method based on detecting a
curated set of protein domains.

Here, we combined the distinct approaches used in previous studies to
create an open-source software package for predicting phage lifestyles
from genome sequence data. `BACPHLIP` (BACterioPHage LIfestyle Predictor)
is a python library with an optional command-line interface that relies
on the HMMER3 software suite [@eddy_accelerated_2011] to identify the
presence of a set of lysogeny-associated protein domains. `BACPHLIP`
assumes that the input genome (nucleotide) sequence is from a fully
complete phage. The phage is initially assumed to be virulent but the
presence and pattern of specific lysogeny-associated protein domains can
override this assumption and result in a temperate classification.

# Development and implementation

We developed `BACPHLIP` by searching the Conserved Domain Database
[@lu_cddsparcle_2020] \(accessed on 03/2020) for protein domains that are
hypothesized to be enriched in temperate phages (*i.e.* mechanistically
involved in lysogeny, see Supplementary Text). We did not include a
broad set of protein domains in our search strategy in order to ensure
interpretability of our model and to limit the possibility of
over-fitting. To determine which of the 371 initial protein domain hits
preferentially associate with temperate phages, we leveraged 1,057
phages with annotated lifestyles collected by
@mavrich_bacteriophage_2017. For each genome sequence, we created a list
of all possible 6-frame translation products $\geq$ 40 amino acids.
Next, we used `HMMER3` to search for the presence of the aforementioned
protein domains, resulting in a vector for each phage describing the
presence (1) or absence (0) of each domain.

At this stage, we randomly split the phage dataset into training and
testing sets (60:40 split, 634 and 423 phages). Using only the training
data, we removed any protein domain that was present in two or fewer
genomes or which was more prevalent in the virulent phage genomes. We
thus established a condensed dataset of 206 putatively useful protein
domains for downstream phage classification. Finally, we fit a Random
Forest classifier to our labeled training data using cross-validation to
tune hyper-parameters (20 separate randomly selected validation sets
drawn from within the training set, see Supplementary Text). The best
performing model from this search, when re-fit to the entire training
set, achieved 99.8% predictive accuracy (633/634 correct predictions) on
the training data.

# Development and implementation

When applied to an independent test set of 423 phages (240 temperate and
183 virulent, drawn randomly from the initial dataset and withheld for
the entirety of model development), BACPHLIP achieved a 98.3%
classification accuracy (415/423 correct predictions). On the same set
of phages, this accuracy exceeded that of both PHACTS and
@mavrich_bacteriophage_2017 (79% and 95.5% accuracy, respectively).
Additional measures of model performance including confusion matrices
for each method are available in the Supplementary Text alongside a
discussion of feature-level contributions (see also Supplementary Table
1 for a complete list of all protein domains).

To ascertain the effect of phylogenetic structure in our dataset (which
could inflate accuracy metrics), we clustered all labeled phages
([@jain_high_2018]) and evaluated accuracy separately on testing set
data for which there were no genomes in the training set with $>$80%
sequence identity across $>$80% of the genome. Although the training set
contained only distant phylogenetic relatives to these 172 genomes,
BACPHLIP still achieved 96.5% accuracy (82% and 93% for PHACTS and
@mavrich_bacteriophage_2017, respectively).

# Discussion

Across all tested datasets, BACPHLIP substantially outperforms existing
methods for classifying phage lifestyles. Our approach relies on access
to training data (phages with known lifestyles) and we encourage users
to recognize the contribution made by @mavrich_bacteriophage_2017. We
emphasize that the existing dataset is made up almost exclusively of
phages from within the *Caudovirales* order and is further biased
towards a small number of hosts (95% infect species within the orders
*Actinobacteria*, *Gammaproteobacteria*, and *Bacilli*). We thus urge
caution when predicting the lifestyle of phages outside of these orders.
We also further reiterate that BACPHLIP was developed for use on
complete phage genomes and performance on fragmented or partially
assembled genomes is likely to be substantially degraded; users are
strongly encouraged to ensure that the starting assumptions are met
prior to running BACPHLIP. We anticipate that the accuracy of BACPHLIP
will increase in future releases as: i) more phylogenetically diverse
phages become available for training the classifier (potentially via
analysis of prophages and/or meta-genomic studies) and ii) discovery and
annotation of conserved protein domains improves to encapsulate new
domains and more phylogenetic diversity amongst the existing protein
domains that BACPHLIP currently relies on.

# Acknowledgements and funding

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References

This work was supported by National Institutes of Health grants F32
GM130113 to A.J.H. and R01 GM088344 to C.O.W.
