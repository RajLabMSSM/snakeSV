#!/bin/bash

snakemake --configfile config/config.yaml --use-conda --conda-create-envs-only --cores 2
