---
title: 'BACPHLIP: Predicting bacteriophage lifestyle from conserved protein domains'
tags:
  - Python
  - bioinformatics
  - genomics
  - bacteriophage
  - microbiology
authors:
  - name: Adam J. Hockenberry^[Corresponding author: adam.hockenberry@utexas.edu]
    orcid: 0000-0003-0872-7098
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Claus O .Wilke
    orcid: 0000-0002-7470-9261
    affiliation: 1
affiliations:
 - name: The University of Texas at Austin
   index: 1
date: 13 August 2017
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
[@paez-espino_uncovering_2016]. 




`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

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
