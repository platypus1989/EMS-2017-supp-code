# Supplementary Code for: Forecasting Urban Household Water Demand with Statistical and Machine Learning Methods Using Large Space-Time Data: A Comparative Study

This repository provides an example data set and code for fitting the models in the 2017 EMS paper. The C++ code provides the functions necessary for training the ST model and using the model to make forecasts. The R scripts simulate some example data and fit and validate the corresponding model. Additional scripts give examples of using the tree-based machine learning methods described in the paper.

## Getting Started

The following instructions will get you a copy of the project up and running on your local machine. You can then change simulation settings or use the model on your own data set.

### Prerequisites

You will need:
* [R](https://cran.r-project.org/), with packages `caret`, `fields`, `mgcv`, `Rcpp`, and `RcppArmadillo` installed
* [RTools](https://cran.r-project.org/bin/windows/Rtools/) (for Windows; other systems may need other compilers, e.g. [here](http://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-os-x/))

### Installing

Once the prerequisites above have been installed, you simply need to copy the files in this repo to your local machine.

### Running the code

The C++ code [`ST-model.cpp`](https://github.com/hrmerrill/EMS-2017-supp-code/blob/master/ST-model.cpp) provides the functions necessary for training the space-time model and using the model to make forecasts. The R script [`ST-model-example-data.R`](https://github.com/hrmerrill/EMS-2017-supp-code/blob/master/ST-model-example-data.R) simulates some example data and fits and validates the ST model. The R scripts [`GAM-model-example-data.R`](https://github.com/hrmerrill/EMS-2017-supp-code/blob/master/GAM-model-example-data.R) and [`machine-learning-example.R`](https://github.com/hrmerrill/EMS-2017-supp-code/blob/master/machine-learning-example.R) simulate some example data and fit and validate the corresponding model (either GAM, or the tree-based machine learning methods described in the paper).

This code should be able to be run as-is, either line-by-line or all at once. For example, running all of [`ST-model-example-data.R`](https://github.com/hrmerrill/EMS-2017-supp-code/blob/master/ST-model-example-data.R) will simulate data and fit and validate the ST model on the simulated data.

## Authors

* **Hunter R. Merrill** - [Personal Site](https://sites.google.com/site/hreidmerrill)
* **Isaac Duerr**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

