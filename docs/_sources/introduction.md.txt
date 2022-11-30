# Introduction
![](D3TaLES_api_logo.png)


API designed for the [D<sup>3</sup>TaLES project](https://d3tales.as.uky.edu/). Includes data processing, database access functions, and property calculators.


## Modules 
This API contains three modules: 
* `D3database`: Interface for interacting with the D<sup>3</sup>TaLES databases. Includes 
queries, data insertion with validation, and a REST API interface.
* `Processors`: Processing code for extracting data from computational and experimental
data files as well as [ExpFlow](https://d3tales.as.uky.edu/expflow) run objects. 
* `Calculators`: Property calculators and plotters for cyclic voltammetry
