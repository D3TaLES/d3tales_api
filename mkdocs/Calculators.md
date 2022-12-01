# Calculators

D<sup>3</sup>TaLES API Calculators allow users to calculate useful experimental and 
computational properties from nested data. Additionally, all calulators contain unit conversion
features. 

Full documentation can be found [here](d3tales_api.Calculators.html).

### Calculator Submodules

The Calculators module contains two primary submodules (`calculators` and `plotters`) and a minor submodules
(`utils`). Most users will only the primary submodules, specifically the `calculators` primary submodule. 

Submodules: 
* `calculators`: This contains many calculator classes. Each calculator class must receive 
`data` and a data `connector` that connects the input data with the connection points. 
Generally, each calculator calculates a specific property with a `calculate` method. 
* `plotters`: The submodule is adapted from the calculators submodule and focuses on producing 
plotting data. Most plotters have a `plot_data` method that produces data for a [plotly](https://plotly.com/)
plot and a `live_plot` method that produces a [Matplotlib](https://matplotlib.org/) plot.  
* `utils`: This contains several useful functions used in the primary submodules, most notable
of which is the `unit_conversion` function. 

## Basic Use



```python

```

The `plotter` submodule has a very similar use to the `calculators` submodule. 

## Considering Units

