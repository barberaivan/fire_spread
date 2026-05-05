# Fire simulation

## Overall structure

This directory contains data and code for a scientific project, that will have a production side. Before keep the repo growing, I want to make a fresh new repo that can easily be expanded in the scientific and production ways.\

- `data` contains raw data, many heavy images
- `files` is made of exports, may work as data, but are files created in this repo
- `dump` has trash, mostly old file that i did not want to lose

Then, remaining folders have mostly code that does things.

3 folders are about fire models-simulation

- `spread` is probably the largest, with all related to fitting the fire spread model. This has all the code for a paper im writing. The spread engine is a C++ R package, FireSpread, but this codes fits the model to data. It generates output that will be constants to use the simulator. So its output is a constant for production, and the result of the science side. It contains the `landscapes_preparation.R`, which runs a loop preparing the landscape that is input for the spread model. This should be a function that can create any landscape, not a loop. It also creates the landscape for the national park (PNNH), which is used in the fire regime simulation.
- `ignition-escape` define and fit the ignition and escape models. The _FWIZ is the most updated; the remaining ones should go away.
- `fire regime simulations` integrates the ignition, escape and spread models, consolidating a fire regime simulator. It has the code that merges models (imports fitted params and all) and here I run many simulations that are scientific side. It also recallibrates some parameters of the spread model.

2 folders related to data:

- `weather` has code that processes weather inputs and exports files. Mostly used for fire regime simulation, but maybe also for spread.
- `flammability indices` has code that fits models that define how raw data will be transformed into fire-models inputs, creating flammability indices. 

The data and files folders have subfolders related to the other main folders, like input-output for-from spread, weather, fire regime, etc.

Files .txt and .md are very useful to understand what happened, what was done, why files are named in certain ways.