# hdnom

[![Build Status](https://travis-ci.org/road2stat/hdnom.svg?branch=master)](https://travis-ci.org/road2stat/hdnom)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/road2stat/hdnom?branch=master&svg=true)](https://ci.appveyor.com/project/road2stat/hdnom)
[![CRAN Version](http://www.r-pkg.org/badges/version/hdnom)](https://cran.r-project.org/package=hdnom)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/hdnom)](http://cranlogs.r-pkg.org/badges/hdnom)

`hdnom` builds nomograms for high-dimensional data with penalized Cox models.

## Paper Citation

Formatted citation for [the preprint](http://dx.doi.org/10.1101/065524):

Nan Xiao, Qing-Song Xu, and Miao-Zhu Li. "hdnom: Building Nomograms for Penalized Cox Models with High-Dimensional Survival Data." bioRxiv (2016): 065524; doi: http://dx.doi.org/10.1101/065524

BibTeX entry:

```
@article {hdnompreprint2016,
	author = {Xiao, Nan and Xu, Qing-Song and Li, Miao-Zhu},
	title = {hdnom: Building Nomograms for Penalized Cox Models with High-Dimensional Survival Data},
	year = {2016},
	doi = {10.1101/065524},
	publisher = {Cold Spring Harbor Labs Journals},
	URL = {http://biorxiv.org/content/early/2016/08/23/065524},
	eprint = {http://biorxiv.org/content/early/2016/08/23/065524.full.pdf},
	journal = {bioRxiv}
}
```

## Installation

To download and install `hdnom` from CRAN:

```r
install.packages("hdnom")
```

Or try the development version on GitHub:

```r
# install.packages("devtools")
devtools::install_github("road2stat/hdnom")
```

To load the package in R, simply use

```r
library("hdnom")
```

and you are all set. See [the vignette](http://hdnom.org/doc/) (can also be opened with `vignette("hdnom")` in R) for a quick-start guide.

## Links

* hdnom project: [http://hdnom.org](http://hdnom.org)
* hdnom web app: [http://hdnom.io](http://hdnom.io)
* hdnom appmaker: [http://hdnom.org/appmaker](http://hdnom.org/appmaker)
* hdnom on CRAN: [https://cran.r-project.org/package=hdnom](https://cran.r-project.org/package=hdnom)
* hdnom on GitHub: [https://github.com/road2stat/hdnom](https://github.com/road2stat/hdnom)
