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
Bacteriophages (viruses that infect bacterial species) are broadly classified
into two distinct lifestyles: temperate phages are capable of a latent phase of
infection within a host cell, whereas virulent phages directly replicate and
lyse host cells upon infection. Determining the lifestyle of a given phage is
critical for understanding its biological significance, but far too many phages
have been discovered in recent years to rely on experimental determination of
this important trait. Here, we present `BACPHLIP`, a computational tool for
predicting bacteriophage lifestyle based solely on genome sequence data.
`BACPHLIP` is a `python` package with an optional command-line interface with
an extensively tuned and tested random forest classifier at its core. On an
independent test set of 423 phage genomes, `BACPHLIP` achieved a classification
accuracy of 98%, greatly exceeding that of the best previously available
software. 

# Statement of need

Bacteriophages play important ecological roles [@paez-espino_uncovering_2016;
@nishimura_environmental_2017; @emerson_host-linked_2018; @daly_viruses_2019;
@gregory_marine_2019], influence both the physiology and evolution of host
species [@touchon_genetic_2016; @carey_phage_2019], and possess a number of
unique traits that are relevant to ongoing biotechnological and medical
applications [@ofir_contemporary_2018; @dedrick_engineered_2019;
@rodriguez-gonzalez_quantitative_2020]. The diversity and availability of phage
genome sequences has expanded in recent years due to the development of
metagenomic sequencing techniques [@deng_viral_2014; @pope_whole_2015;
@roux_viral_2015; @tisza_discovery_2020]. A particularly important aspect of
bacteriophage biology is whether individual species have a temperate (lysogenic) or
virulent (lytic) lifestyle [@bobay_adaptation_2013], a classification which
refers to the ability of some phage species to integrate into host bacterial cell
genomes where they may lay dormant for extended periods of time.

The gold-standard for phage lifestyle classification relies on experimental
determination, but these approaches are impractical for the entirety of newly
discovered phage genomes (many of whom are known *only* by their genomes).
Using computational methods, @mcnair_phacts_2012 developed a random forest
classifier (`PHACTS`) to predict lifestyle based off of sequence similarity to
a set of query proteins. More recently, @mavrich_bacteriophage_2017 described a
computational classification method based on detecting a curated set of protein
domains but did not make this software available.

Here, we combined the distinct approaches used in previous studies to create an
open-source software package for predicting phage lifestyles from genome
sequence input (`fasta` formatted). `BACPHLIP` (BACterioPHage LIfestyle
Predictor) is a `python` library with an optional command-line interface that
relies on the `HMMER3` software suite [@eddy_accelerated_2011] to identify the
presence of a set of lysogeny-associated protein domains within the genome
(methodically selected from the `Conserved Domain Database`
[@lu_cddsparcle_2020]). Using this presence/absence information, the core of
`BACPHLIP` is a random forest classifier that reports on the probability that
the input phage genome is either temperate or virulent. This approach achieved
a classification accuracy of 98.3% (415/423 correct predictions) on a fully
independent set of testing data, greatly exceeding that of the existing
`PHACTS` software (79%). 

# Modeling assumptions

We emphasize that the `BACPHLIP` model relies on an existing labeled dataset that
is made up almost exclusively of phages from within the *Caudovirales* order
and is further biased towards a small number of hosts (95% infect species
within the orders *Actinobacteria*, *Gammaproteobacteria*, and *Bacilli*). We
thus urge caution when predicting the lifestyle of phages outside of these
phylogenetic orders. We also note that `BACPHLIP` was developed
for use on complete phage genomes and performance on fragmented or partially
assembled genomes is likely to be degraded; users are strongly encouraged to
ensure that the starting assumptions (re-itereated extensively in the package
documentation) are met prior to using `BACPHLIP`. We anticipate that the
accuracy of `BACPHLIP` will increase in future releases as: i) more
phylogenetically diverse phages become available for training the classifier
(potentially via analysis of prophages and/or meta-genomic studies) and ii)
discovery and annotation of conserved protein domains improves to encapsulate
new domains and more phylogenetic diversity amongst the existing protein
domains that `BACPHLIP` currently relies on.

# Acknowledgements

This work was supported by National Institutes of Health grants F32
GM130113 to A.J.H. and R01 GM088344 to C.O.W.

# References

