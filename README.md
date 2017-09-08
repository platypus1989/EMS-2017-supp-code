# Supplementary Code for: Forecasting Urban Water Demand with Statistical and Machine Learning Methods Using Large Space-Time Data

This repository provides an example data set and code for fitting the spatio-temporal model in the 2017 EMS paper. The C++ code provides the functions necessary for training the model and using the model to make forecasts. The R code simulates some example data, fits a linear regression model and the spatio-temporal model, and compares forecasts made by each model.

## Getting Started

The following instructions will get you a copy of the project up and running on your local machine. You can then change simulation settings or use the model on your own data set.

### Prerequisites

You will need:
* [R](https://cran.r-project.org/), with packages `fields`, `Rcpp`, and `RcppArmadillo` installed
* [RTools](https://cran.r-project.org/bin/windows/Rtools/) (for Windows; other systems may need other compilers, e.g. [here](http://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-os-x/))

### Installing

Once the prerequisites above have been installed, you simply need to copy the two files in this repo ([ST-model.cpp](ST-model.cpp) and [ST-model-example-data.R](ST-model-example-data.R)) to your local machine.

### Running the code

This code should be able to be run as-is, either line-by-line or all at once. 

## Authors

* **Hunter R. Merrill** - [Personal Site](https://sites.google.com/site/hreidmerrill)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

