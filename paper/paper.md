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
on the HMMER3 software suite ([@eddy_accelerated_2011]) to identify the
presence of a set of lysogeny-associated protein domains. `BACPHLIP`
assumes that the input genome (nucleotide) sequence is from a fully
complete phage. The phage is initially assumed to be virulent but the
presence and pattern of specific lysogeny-associated protein domains can
override this assumption and result in a temperate classification.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
