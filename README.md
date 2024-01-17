# GLMsim
GLMsim is a single cell simulator that can simultaneously capture the library size, biology and unwanted variation and their associations via a generalized linear model, and to simulate data resembling the original experimental data in these respects. GLMsim is capable of quantitatively benchmarking different single cell integration methods, and assessing their abilities to retain biology and remove library size and batch effects.
## Installation
`devtools::install_github("jiananwehi/GLMsim")`
## Preliminary setup
The following packages are prerequisites to run GLMsim. <br />
`library(GLMsim)`
`library(MASS)` 
`library(parallel)`
