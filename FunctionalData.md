---
name: FunctionalData
topic: Functional Data Analysis
maintainer: Fabian Scheipl, Eleonora Arnone, Giles Hooker, J. Derek Tucker, Julia Wrobel
email: fabian.scheipl@stat.uni-muenchen.de
version: 2024-06-17
source: https://github.com/cran-task-views/FunctionalData/
---

Functional data analysis (FDA) deals with data that ["provides
information about curves, surfaces or anything else varying over a
continuum."](https://en.wikipedia.org/wiki/Functional_data_analysis)
This task view tries to provide an overview of available packages in this developing
field.  
In practice, there is substantial overlap between time series data and functional data,
so many packages listed under the `r view("TimeSeries")` task view
could be useful for functional data as well and vice versa.

### General functional data analysis

The packages listed below provide "infrastructure" for representing and handling
function-valued data and/or implement many widely applicable functional
data methods:

-   `r pkg("fda", priority = "core")` provides object-types for functional data 
    with corresponding functions for smoothing, plotting and simple regression models,
    c.f. Ramsay et al. (2009, `r doi("10.1007/978-0-387-98185-7")`).
-   `r pkg("fdapace", priority = "core")` provides functional
    principal component based methods for sparsely or densely sampled
    random trajectories and time courses for functional regression and
    correlation, for longitudinal data analysis, the analysis of
    stochastic processes from samples of realized trajectories, and for
    the analysis of underlying dynamics.
-   `r pkg("fdasrvf")` performs alignment, PCA, and
    regression of multidimensional or unidimensional functions using the
    square-root velocity framework (Srivastava et al., 2011). This
    framework allows for elastic analysis of functional data through
    phase and amplitude separation.
-   `r pkg("fda.usc", priority = "core")` provides routines
    for exploratory and descriptive analysis of functional data such as
    depth measurements, outlier detection, as well as unsupervised and
    supervised classification, (univariate, nonparametric) regression
    models with a functional covariate and functional analysis of
    variance.
-   `r pkg("fds", priority = "core")` contains 19 data sets
    with functional data.
-   `r pkg("funData")` provides S4 classes for univariate
    and multivariate functional and image data and utility functions.
-   `r pkg("tf", priority = "core")` provides S3 vector data types for functional data with arithmetic and summary methods, derivation, integration and smoothing, plotting, data import and export, and data wrangling (subsetting, sub-assigning, zooming, extracting functional features) that allow including such vectors in data frames for joint analysis of functional and scalar variables.

### Inference methods for functional data

**General purpose regression and classification:**

-   `r pkg("FDboost", priority = "core")` implements flexible
    additive regression models and variable selection for
    scalar-on-function, function-on-scalar and function-on-function
    regression models that are fitted by a component-wise gradient
    boosting algorithm.
-   `r pkg("mlr3fda")` extends ['mlr3'](https://mlr3.mlr-org.com/) for functional data analysis, with support for irregular and regular functional data as defined in ['tf'](https://tidyfun.github.io/tf/), providing 'PipeOps' for preprocessing functional columns and for extracting scalar features to be used in standard machine learning algorithms.
-   `r pkg("refund", priority = "core")` provides
    spline-based methods for penalized function-on-scalar,
    scalar-on-function, and function-on-function regression as well as
    methods for functional PCA. Some of the functions are also applicable to
    image data.
-   `r pkg("dbstats")` provides prediction methods where
    explanatory information is coded as a matrix of distances between
    individuals. It includes distance based versions of
    `lm` and `glm`, as well as nonparametric versions of both, based on local estimation.
    To apply these methods to functional data it is sufficient to
    calculate a distance matrix between the observed functional data.



**Specialized models and applications:**

-   `r pkg("cfda")` implements analysis of categorical functional data, by computing an encoding (a real functional variable) for each state.
-   `r pkg("denseFLMM")` and `r pkg("sparseFLMM")` estimate linear mixed
    models for densely and sparsely sampled functional responses, respectively, based on
    functional principal component analysis.
- `r pkg("fastFMM")` implements a fast univariate inference approach for fitting functional mixed models.   
-   `r pkg("fdaconcur")` implements concurrent and history index models for both dense and sparsely observed functional data
-   `r pkg("fdaPDE")` implements statistical analysis of functional and spatial 
   data over multidimensional complex domains, based on regression models with 
   partial differential regularization, discretized through the finite element method.
-   `r pkg("GPFDA")` uses functional regression as the mean
    structure and Gaussian processes as the covariance structure.
-   `r pkg("growfunctions")` estimates a collection of
    time-indexed functions under either of Gaussian process (GP) or
    intrinsic Gaussian Markov random field (iGMRF) prior formulations
    where a Dirichlet process mixture allows sub-groupings of the
    functions to share the same covariance or precision parameters. 
-   `r pkg("multifamm")` implements multivariate functional additive mixed models
    (multiFAMM) based on univariate sparse functional regression models and 
    multivariate functional principal component analysis, c.f. 
    Volkmann et al. (2012, `r doi("10.48550/arXiv.2103.06606")`).
-   `r pkg("robflreg")` implements robust methods for functional linear regression
    (scalar-on-function, function-on-function) based on robust partial
    least squares decompositions.
- `r pkg("slasso")` implements a smooth LASSO estimator for the function-on-function linear regression model, 
   c.f. Centofanti et al. (2020, `r doi("10.48550/arXiv.2007.00529")`).     
-   `r pkg("splinetree")` implements regression trees and
    random forests for longitudinal or functional data using a spline
    projection method.
-   `r pkg("adass")` implements the adaptive smoothing spline estimator for the function-on-function linear regression model described, c.f. Centofanti et al. (2023, `r doi("10.1007/s00180-022-01223-6")`).
  
**Statistical tests:**

-   `r pkg("fdANOVA")` implements analysis of variance testing procedures for univariate and multivariate functional data.
-   `r pkg("fdatest")` provides an implementation of the Interval Testing Procedure for functional data in different frameworks (i.e., one or two-population frameworks, functional linear models) by means of different basis expansions (i.e., B-spline, Fourier, and phase-amplitude Fourier).
-   `r pkg("gmtFD")` implements general multiple tests for univariate and multivariate functional data.
-   `r pkg("multiFANOVA")`  implements multiple contrast tests for functional data for one-way and multi-way analysis of variance for functional data.
-  `r pkg("rmfanova")`  implements statistical tests for functional repeated measures analysis testing for equality or pairwise comparisons (post-hoc analysis) of mean functions in repeated experiments.
- `r pkg("rofanova")`  implements robust functional analysis of variance, c.f. Centofanti et al. (2021, `r doi("10.48550/arXiv.2112.10643")`).
-   `r pkg("SCBmeanfd")` provides methods for estimating and inferring the mean of functional data. The methods include simultaneous confidence bands, local polynomial fitting, bandwidth selection by plug-in and cross-validation, goodness-of-fit tests for parametric models, equality tests for two-sample problems, and plotting functions.

### Clustering functional data

-   `r pkg("fdacluster")` implements k-means, hierarchical agglomerative and DBSCAN clustering for functional data which allows for jointly aligning and clustering curves.
-   `r pkg("funFEM")` 's algorithm (Bouveyron et al., 2014)
    allows to cluster functional data by modeling the curves within a
    common and discriminative functional subspace.


### Registering and aligning functional data

-   `r pkg("elasdics")` provides functions to align 2D curves and to compute mean curves based on the elastic distance defined in the square-root-velocity 
   framework, c.f. Steyer et al. (2021, `r doi("10.48550/arXiv.2104.11039")`).
-   `r pkg("elastes")` provides functions for the computation of functional elastic shape means over sets of open planar curves,  particularly for settings where these curves are only sparsely and irregularly observed, c.f. Stöcker et al. (2022, `r doi("10.48550/arXiv.2203.10522")`). 
-   `r pkg("fdasrvf")` performs alignment, PCA, and
    regression of multidimensional or unidimensional functions using the
    square-root velocity framework (Srivastava et al., 2011). This
    framework allows for elastic analysis of functional data through
    phase and amplitude separation.
-   `r pkg("registr")` provides registration for (incomplete) non-Gaussian
    functional data,  c.f Wrobel et al. (2019, `r doi("10.1111/biom.12963")`),
    Wrobel and Bauer (2021, `r doi("10.21105/joss.02964")`).
-   `r pkg("warpMix")` implements warping (alignment) for
    functional data using B-spline based mixed effects models.

### Time series of functional data

-   `r pkg("ftsa", priority = "core")` provides functions for
    visualizing, modeling, forecasting and hypothesis testing of
    functional time series.
-   `r pkg("ftsspec")` provides functions for estimating the
    spectral density operator of functional time series (FTS) and
    comparing the spectral density operator of two functional time
    series, in a way that allows detection of differences of the
    spectral density operator in frequencies and along the curve length.
-   `r pkg("freqdom")` provides frequency domain methods for
    multivariate and functional time series and implements dynamic
    functional principal components and functional regression in the
    presence of temporal dependence.
-   `r pkg("fdaACF")` contains functions to quantify the
    serial correlation across lags of a given functional time series, see also 
    [github.com/GMestreM/fdaACF](https://github.com/GMestreM/fdaACF) .

### Visualization and Exploratory Data Analysis

-   `r pkg("rainbow")` contains functional data sets and functions for 
    functional data display, exploratory analysis and outlier detection.
-   `r pkg("refund.shiny")` provides a function `plot_shiny` that uses S3 and generates a Shiny app of interactive visualizations for functional data analyses produce using the `refund` package, described in c.f. Wrobel et al. (2016, `r doi("10.1002/sta4.109")`)
-   `r pkg("fdaoutlier")` provides a collection of functions
    for functional data outlier detection. Methods implemented include
    directional outlyingness, MS-plot, total variation depth, and
    sequential transformations among others.

### Other

- `r pkg("conformalInference.fd")` computes full conformal, split conformal and multi split conformal prediction regions when the response has functional nature. 
-   `r pkg("ddalpha")` implements depth-based classification
    and calculation of data depth, also for functional data.
-   `r pkg("face")` implements Fast Covariance Estimation
    for Sparse Functional Data, c.f. Xiao et al. (2018, `r doi("10.1007/s11222-017-9744-8")`).
-   `r pkg("fdadensity")` implements Petersen and Müller (2016, `r doi("10.1214/15-AOS1363")`) for the
    analysis of samples of density functions via specialized Functional
    Principal Components Analysis.
-    `r pkg("fdaPOIFD")` implements partially observed integrated functional depth with applications to visualization, outlier detection and classification, c.f.  Elías et al. (2022, `r doi("10.1080/10618600.2022.2070171")`).
-   `r pkg("frechet")` implements Fréchet regression for for non-Euclidean responses, 
     e.g. distributions in L^2-Wasserstein space or covariance matrices, 
     c.f. Petersen & Müller (2019, `r doi("10.1214/17-AOS1624")`).
- `r pkg("funcharts")` provides functional control charts for statistical process monitoring of functional data, c.f. Capezza et al (2023, `r doi("10.1080/00224065.2023.2219012")`)     
-   `r pkg("mfaces")` implements multivariate functional
    principal component analysis via fast covariance estimation for
    multivariate sparse functional data or longitudinal data,  c.f Li,
    Xiao, and Luo (2020, `r doi("10.1002/sta4.245")`) 
-   `r pkg("MFPCA")` calculates multivariate FPCA for "multimodal" data 
    observed on domains with different dimensionalities, c.f. Happ and Greven (2018, `r doi("10.1080/01621459.2016.1273115")`).



If you want to suggest additions or improvements to this list,
please submit an issue or pull request in the 
[GitHub repository](https://github.com/cran-task-views/FunctionalData/).



### Links

-   [Website of the canonical FDA book by Ramsay and Silverman](http://www.psych.mcgill.ca/misc/fda/)
-   [PACE: collection of MATLAB scripts from UC Davis](http://www.stat.ucdavis.edu/PACE/)
-   [WFMM: powerful software for Bayesian wavelet-based functional mixed models (C++/Matlab)](https://biostatistics.mdanderson.org/softwaredownload/SingleSoftware.aspx?Software_Id=70)
-   [scikit-FDA: comprehensive Python package for FDA](https://fda.readthedocs.io)
-   [fdasrsf: comprehensive Python package for SRVF FDA](https://fdasrsf-python.readthedocs.io/)
