# Interactive-metal-mixture-toxicity-as-emergent-property

## Context & Purpose

The code provided in this repository is developed to model and predict metal mixture toxicity to *Daphnia magna* populations in laboratory microcosms.
See article published in [ET&C](https://setac.onlinelibrary.wiley.com/doi/full/10.1002/etc.5176). 

## Usage

The function `Run` can be used to execute a simulation of *D. magna* population dynamics. It takes two dictionaries,
`global_dict` and `fleas_dict` as arguments. These dictionaries contain the global and individual-level simulation parameters. Optional keyword arguments are Arrays of exposure concentrations for each metal, which can be used to simulate time-varying exposure. <br>
`SingleStressorTest` and `MixtureTest` are functions to simulate toxicity tests. Under the hood, these are just repeated calls to `Run`. <br>
Concrete examples are given in the [figshare repository](https://figshare.com/articles/journal_contribution/Interactive_metal_mixture_toxicity_to_D_magna_populations/14161283) associated with the ET&C article.
