# Description

This directory contains the files required for building and running the transcription submodule.

# Process

The main process in transcription is a Markov chain describing the states of the polymerase. These states can be
* Free (unbound, floating in cytosol)
* Non-specifically bound to the chromosome (exact locus is not recorded)
* Specifically bound to the chromosome at a certain transcription unit
* Actively transcribing, implies that the polymerase is actively elongating the nucleotide chain

Transcription factors (sigma factor) are required for the specifically bound -> active transition.

# Useage

From the root of the repository, the Python script can be run as follows:
```
python src/transcription/map_tx_unit_to_model.py --tx-units data/s3k-transcription-units.csv --genes data/genes.csv --locus --tu-cutoff 1 --bp-cutoff 2
```

This generase a temporary file (`/tmp/tx.sbml`) which contains the SBML model of the transcription process. This model can be run using Tellurium 1.1.8 using the file `run_model.py`. Simulations were carried out on Mac OS X 10.10.