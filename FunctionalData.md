---
name: FunctionalData
topic: Functional Data Analysis
maintainer: Fabian Scheipl
email: fabian.scheipl@stat.uni-muenchen.de
version: 2022-03-07
source: https://github.com/cran-task-views/FunctionalData/
---

Functional data analysis (FDA) deals with data that ["provides
information about curves, surfaces or anything else varying over a
continuum."](https://en.wikipedia.org/wiki/Functional_data_analysis)
This task view catalogues available packages in this rapidly developing
field.

### General functional data analysis

-   `r pkg("fda", priority = "core")` provides functions to
    enable all aspects of functional data analysis: It includes
    object-types for functional data with corresponding functions for
    smoothing, plotting and regression models. The package includes data
    sets and script files for working examples from the book: Ramsay, J.
    O., Hooker, Giles, and Graves, Spencer (2009) "Data Analysis with R
    and Matlab" (Springer).
-   `r pkg("fdasrvf", priority = "core")` performs alignment,
    PCA, and regression of multidimensional or unidimensional functions
    using the square-root velocity framework (Srivastava et al., 2011).
    This framework allows for elastic analysis of functional data
    through phase and amplitude separation.
-   `r pkg("fdapace", priority = "core")` provides functional
    principal component based methods for sparsely or densely sampled
    random trajectories and time courses for functional regression and
    correlation, for longitudinal data analysis, the analysis of
    stochastic processes from samples of realized trajectories, and for
    the analysis of underlying dynamics.
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
-   `r pkg("rainbow")` contains functions and data sets for
    functional data display, exploratory analysis and outlier detection.
-   `r pkg("fdaoutlier")` provides a collection of functions
    for functional data outlier detection. Methods implemented include
    directional outlyingness, MS-plot, total variation depth, and
    sequential transformations among others.

### Regression and classification for functional data

-   `r pkg("dbstats")` provides prediction methods where
    explanatory information is coded as a matrix of distances between
    individuals. It includes distance based versions of
    `             lm           ` and `             glm,           ` as
    well as nonparametric versions of both, based on local estimation.
    To apply these methods to functional data it is sufficient to
    calculate a distance matrix between the observed functional data.
-   `r pkg("denseFLMM")` and
    `r pkg("sparseFLMM")` estimate functional linear mixed
    models for densely and sparsely sampled data, respectively, based on
    functional principal component analysis.
-   `r pkg("fdANOVA")` implements analysis of variance
    testing procedures for univariate and multivariate functional data
-   `r pkg("FDboost", priority = "core")` implements flexible
    additive regression models and variable selection for
    scalar-on-function, function-on-scalar and function-on-function
    regression models that are fitted by a component-wise gradient
    boosting algorithm.
-   `r pkg("flars")` implements variable selection for the
    functional linear regression with scalar response variable and mixed
    scalar/functional predictors based on the least angle regression
    approach.
-   `r pkg("GPFDA")` uses functional regression as the mean
    structure and Gaussian processes as the covariance structure.
-   `r pkg("growfunctions")` estimates a collection of
    time-indexed functions under either of Gaussian process (GP) or
    intrinsic Gaussian Markov random field (iGMRF) prior formulations
    where a Dirichlet process mixture allows sub-groupings of the
    functions to share the same covariance or precision parameters. The
    GP and iGMRF formulations both support any number of additive
    covariance or precision terms, respectively, expressing either or
    both of multiple trend and seasonality.
-   `r pkg("refund", priority = "core")` provides
    spline-based methods for roughness penalized function-on-scalar,
    scalar-on-function, and function-on-function regression as well as
    methods for functional PCA. Some of the functions are applicable to
    image data.
-   `r pkg("splinetree")` implements regression trees and
    random forests for longitudinal or functional data using a spline
    projection method.

### Clustering functional data

-   `r pkg("funFEM")` 's algorithm (Bouveyron et al., 2014)
    allows to cluster functional data by modeling the curves within a
    common and discriminative functional subspace.
-   `r pkg("funHDDC")` provides the funHDDC algorithm
    (Bouveyron & Jacques, 2011) which allows to cluster functional data
    by modeling each group within a specific functional subspace.
-   `r pkg("funLBM")` implements model-based co-clustering
    of functional data, i.e., simultaneously clustering the rows and the
    columns of a data matrix where each entry of the matrix is a
    function or a time series.
-   `r pkg("fdakma")` performs clustering and alignment of a
    multidimensional or unidimensional functional dataset by means of
    k-mean alignment.

### Registering and aligning functional data

-   `r pkg("fdasrvf")` performs alignment, PCA, and
    regression of multidimensional or unidimensional functions using the
    square-root velocity framework (Srivastava et al., 2011). This
    framework allows for elastic analysis of functional data through
    phase and amplitude separation.
-   `r pkg("fdakma")` performs clustering and alignment of a
    multidimensional or unidimensional functional dataset by means of
    k-mean alignment.
-   `r pkg("registr")` provides registration for
    (incomplete) non-Gaussian functional data, c.f Wrobel et al. (2019)
    [doi: 10.1111/biom.12963](https://doi.org/10.1111/biom.12963) and
    Wrobel and Bauer (2021) [doi:
    10.21105/joss.02964](https://doi.org/10.21105/joss.02964) .
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
-   `r pkg("freqdom.fda")` provides a wrapper for
    functionality of `r pkg("freqdom")` for objects from
    `r pkg("fda")`
-   `r pkg("pcdpca")` extends multivariate dynamic principal
    components to periodically correlated multivariate and functional
    time series.
-   `r pkg("fdaACF")` contains functions to quantify the
    serial correlation across lags of a given functional time series
    using an autocorrelation function for functional time series. The
    autocorrelation function is based on the L2 norm of the lagged
    covariance operators of the series. Functions are available for
    estimating the distribution of the autocorrelation function under
    the assumption of strong functional white noise. A brief
    illustration of the functionality of the proposed functions can be
    seen at
    [github.com/GMestreM/fdaACF](https://github.com/GMestreM/fdaACF) .

### Other

-   `r pkg("covsep")` provides functions for testing if the
    covariance structure of 2-dimensional data is separable.
-   `r pkg("ddalpha")` implements depth-based classification
    and calculation of data depth, also for functional data.
-   `r pkg("fdadensity")` implements Petersen and
    Mueller (2016) ( [doi:
    10.1214/15-AOS1363](https://doi.org/10.1214/15-AOS1363) ) for the
    analysis of samples of density functions via specialized Functional
    Principal Components Analysis.
-   `r pkg("face")` implements Fast Covariance Estimation
    for Sparse Functional Data paper (c.f. Statistics and Computing,
    [doi:
    10.1007/s11222-017-9744-8](https://doi.org/10.1007/s11222-017-9744-8)
    )
-   `r pkg("fdatest")` provides an implementation of the
    Interval Testing Procedure for functional data in different
    frameworks (i.e., one or two-population frameworks, functional
    linear models) by means of different basis expansions (i.e.,
    B-spline, Fourier, and phase-amplitude Fourier).
-   `r pkg("fdcov")` provides a variety of tools for the
    analysis of covariance operators.
-   `r pkg("geofd")` provides Kriging based methods for
    predicting functional data (curves) with spatial dependence.
-   `r pkg("mfaces")` implements multivariate functional
    principal component analysis via fast covariance estimation for
    multivariate sparse functional data or longitudinal data (c.f Li,
    Xiao, and Luo (2020) [doi:
    10.1002/sta4.245](https://doi.org/10.1002/sta4.245) ).
-   `r pkg("MFPCA")` calculates multivariate FPCA for data
    observed on different dimensional domains (c.f. Happ and Greven
    (2018) [doi:
    10.1080/01621459.2016.1273115](https://doi.org/10.1080/01621459.2016.1273115)
    ).
-   `r pkg("SCBmeanfd")` provides methods for estimating and
    inferring the mean of functional data. The methods include
    simultaneous confidence bands, local polynomial fitting, bandwidth
    selection by plug-in and cross-validation, goodness-of-fit tests for
    parametric models, equality tests for two-sample problems, and
    plotting functions.
-   `r pkg("switchnpreg")` provides functions for estimating
    the parameters from the latent state process and the functions
    corresponding to the J states as proposed by De Souza and Heckman
    (2013).

Please e-mail the maintainer with suggestions, additions, and
improvements or submit an issue or pull request in the GitHub repository
linked above.


### Links
-   [Website of the canonical FDA book by Ramsay and Silverman](http://www.psych.mcgill.ca/misc/fda/)
-   [PACE: collection of MATLAB scripts from UC Davis](http://www.stat.ucdavis.edu/PACE/)
-   [WFMM: powerful software for Bayesian wavelet-based functional mixed models (C++/Matlab)](https://biostatistics.mdanderson.org/softwaredownload/SingleSoftware.aspx?Software_Id=70)
-   [scikit-FDA: comprehensive Python package for FDA](https://fda.readthedocs.io)
