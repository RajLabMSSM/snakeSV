# Tutorial for installing snakeSV

## Preparing the environment

Let's install `Miniconda` to manage easily install tools in your virtual machine.

`Miniconda` is a free minimal installer for **conda** that includes only `conda`, Python and a small number of other useful packages. `Miniconda` allows you to create a minimal self contained Python installation, and then use the `conda` command to install additional packages.

Run the following commands to install `Miniconda`:
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```

Reload your shell to enable the `conda` command
```bash
bash
```

Add `bioconda` and `conda-forge` to your channels
```bash 
conda config --add channels conda-forge
conda config --add channels bioconda
```

Install `mamba` to your `base` environment
```bash
conda install -y mamba
```

Now create an environment named `snakesv_env` and install `snakeSV` in it.
```bash
mamba create -n snakesv_env snakesv
```

Activate the environment
```bash
conda activate snakesv_env
```

Finally another run the command to test the installation 
```bash
snakeSV --test_run
```

Results can be found in the `results_snakesv` folder.
