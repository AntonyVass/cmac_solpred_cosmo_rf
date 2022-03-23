# cmac_solpred_cosmo_rf
Code used in the generation of results for "A Unified AI Framework for Solubility Prediction Across Organic Solvents"

## Disclaimer
All code is provided AS IS without warranty of any kind. It was written by non-experts in order to perform a task once. Please exercise caution and initiative in deploying this code, and forgive the rampant bad practices visible throughout (hard-coding, copy/pasting sections...)

## Usage
All scripts require the user to set I/O specific to their environment.  We have moved most of these lines to the top in order to ease editing.  If working with the provided dataset (recommended in order to familiarise yourself with our scrappy code), little further editing is needed.

Execution time for single RF models is approx. 1 min on a modern desktop machine.  The drip-feeding loop will vary greatly depending on the number of train/test combinations specified.  The full set of results generated for the manuscript took approx. 190 hours.
