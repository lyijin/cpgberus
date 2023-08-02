# Processing ONT data for this project #

The four- or five-base ONT data (FAST5/FASTQ) produced for this project cannot be uploaded due to potential identifiability concerns. The BED files for the per-position methylation levels (near negligible risk of re-identification) are available at TODO:CSIRO DAP.

## Introduction ##

The cool/good thing about raw ONT data is that you get separate folders for FAST5 (raw raw data) and basecalled FASTQ data, binned into "PASS" or "FAIL" for each type of data.

As relative newcomers to this type of data, we did not second-guess ONT's calls, and used the FAST5 PASS and FASTQ PASS contents for methylation calling and coverage calculations respectively.

The coverage calling from FASTQ is relatively straightforward: use FASTQs and map them against reference sequences. Exact commands are in the `15_ont_minimap2_coverage/` folder.

The direct methylation calling from FAST5 is... not so easy. Read on, but remember, here be dragons. Guidance likely to be out-of-date really quickly due to the speed of ONT updating their bioinformatics tools, but well at least we've done our part in easing replicability of science.

If this is too much of a headache, the scripts in `16_loci_specific_three_way/` will work with the intermediate files we provide in that folder.

## Setting up `megalodon` on a `slurm` cluster ##

We ran `megalodon` on the CSIRO `slurm` cluster.

`megalodon` is an ONT-authored tool that calls modified bases from ONT reads. In this example, we wanted to call 5mC bases from reads produced from a GridION ( = MinION) machine.

`megalodon` depends on `guppy` and also on `remora`, all from ONT.

All three tools are under active development, so this pipeline (last updated Feb 2022) might not be best practice in a year or two.

## Setting up `guppy` ##

Download URL can be obtained from https://community.nanoporetech.com/downloads, but requires login.

However,  once the URL is known, no login credentials to download; i.e. download can be done on the cluster directly.

```shell
cd ~
mkdir tools
cd tools/
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_5.0.14_linux64.tar.gz
tar zxvf ont-guppy_5.0.14_linux64.tar.gz
```

## Setting up `megalodon` ##

From https://github.com/nanoporetech/megalodon, plus modifications.

```shell
module load miniconda3                # if the module isn't on slurm, you'll need to install miniconda3
conda create -n megalodon python=3.8  # atm, ont-pyguppy-client-api needs python < 3.9
conda init bash                       # run this if prompted by conda; remember to restart your shell
conda activate megalodon
pip install megalodon

megalodon --help                      # should produce help screen, not "command not found"
megalodon --version                   # the version on pip might not be new enough to support all the
                                      # fancy shmancy stuff

# to get the absolutely most advanced version possible
git clone https://github.com/nanoporetech/megalodon
pip install -e megalodon/
```

## Setting up `remora` ##

`remora` is in the process of superseding `rerio`, best to get used to using this now and ignore `rerio`.

Docs are at https://github.com/nanoporetech/remora, installation is much easier than `rerio`.

```shell
pip install ont-remora
```

## (deprecated) Setting up `rerio` ##

From asking ONT and reading docs, `rerio` was trained on a mix of WGBS and synthetic datasets (PCR-ed DNA for unmethylated DNA, M.SssI-treated for methylated DNA). `remora` was however trained purely on synthetic datasets.

From https://github.com/nanoporetech/rerio, plus modifications.

```shell
git clone https://github.com/nanoporetech/rerio
# models are poretype-specific and modified-base-specific. check `rerio` github for a full list of models; alter to taste
rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases_5mC_CpG_v001
cp ont-guppy/data/barcoding/* rerio/basecall_models/barcoding/  # enable barcoding support
```

# SLURM script used in this work #

`megalodon` DOES NOT have multiplexing capabilities yet--it won't automatically split results by barcode; this has to be done manually. \
https://github.com/nanoporetech/megalodon/issues/43. 

SLURM script that splits reads by barcode, for input data structured as

```shell
$ ls -l */

00_raw_reads/:
total 4.0K
drwxr-xr-x 7 lie128 hpc-users 4.0K Feb 17 12:15 20210914_0348_X4_FAQ88026_aa298ba1/

03_remora_vs_45s/:
total 4.0K
drwxr-xr-x 5 lie128 hpc-users 4.0K Feb 17 12:28 ont-guppy/

refs/:
total 280K
-rw-r--r-- 1 lie128 hpc-users 278K Jun 17  2021 homo_sapiens.KY962518.last_11kb_first.ont.mmi
```

and

```shell
$ ls -l 00_raw_data/20210914_0348_X4_FAQ88026_aa298ba1/fast5_pass/

total 52K
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode13/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode14/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode15/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode16/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode17/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode18/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode19/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode20/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode21/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode22/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode23/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 barcode24/
drwxr-xr-x 2 lie128 hpc-users 4.0K Sep 16 12:22 unclassified/
```

This SLURM script only calls 5mC in the CpG context, as we weren't interested in non-CpG methylation nor 5hmC bases.

```shell
#!/bin/bash
#SBATCH --job-name=test_meg_gpu
#SBATCH --time=1:59:00
#SBATCH --ntasks-per-node=14
#SBATCH --gres=gpu:1
#SBATCH --mem=16g

# application specific commands
#module load miniconda3
#conda activate megalodon

echo assigned CUDA device: ${CUDA_VISIBLE_DEVICES}

for a in ../00_raw_reads/20210914_0348_X4_FAQ88026_aa298ba1/fast5_pass/*
do
  b=`basename ${a}`
  megalodon ${a} \
    --outputs basecalls mappings mod_mappings mods \
    --output-directory ${b} \
    --reference ../refs/homo_sapiens.KY962518.last_11kb_first.ont.mmi \
    --guppy-server-path ./ont-guppy/bin/guppy_basecall_server \
    --guppy-config dna_r9.4.1_450bps_fast.cfg \  # this file depends on what kit/pore type was used in expt
    --guppy-timeout 600 \
    --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5mc CG 0 \  # this file depends on what kit/pore type was used in expt
    --suppress-progress-bars \
    --devices "${CUDA_VISIBLE_DEVICES}" --processes 14 --overwrite
done
```

Note: `*.mmi` files are minimap2 index files, generated with commands like \
`minimap2 -x map-ont -d homo_sapiens.KY962518.last_11kb_first.ont.mmi homo_sapiens.KY962518.last_11kb_first.fa`

Note 2: GPU `megalodon` was ~100x faster than CPU `megalodon`. Highly recommend getting GPU to work, something that took 3 days to run on CPU took ~40 mins to complete.

# Results #

The resulting files were all named "modified_bases.5mC.bed" files located in folders with the respective barcodes. These files are in bedMethyl format. \
https://www.encodeproject.org/data-standards/wgbs/

TL;DR coverage is in penultimate column, methyl % is in final column.

We renamed these files (with sample names e.g., WR025V1O) so that they could all be placed in the same folder for plotting purposes.

# Postscript #

## Potholes faced ##

We encountered weird compilation errors when the default compiler is `icc` (Intel's compiler). You can find out what the default is by running

```shell
module list
Currently Loaded Modulefiles:
  1) SC   2) slurm/current(default)   3) cuda-driver/current   4) intel-cc/16.0.1.150   5) intel-fc/16.0.1.150
```

If you see `intel-cc` and `intel-fc`, then yes, you're on icc.

Change this to `gcc` to avoid the potholes we encountered when we initially used `icc` (see "Troubleshooting corner" below)

```shell
module unload intel-cc intel-fc

module avail -l gcc
- Package/Alias -----------------------.- Versions --------.- Last mod. -------
/apps/modules/modulefiles:
gcc/9.3.0                               default             2020/11/12 12:35:51
gcc/10.3.0                                                  2020/11/12 13:47:38
gcc/11.1.0                                                  2021/04/29 22:47:27

module load gcc/11.1.0    # i belong to the cult of the new
```

Also, we like to collate downloaded tools in `~/tools`. It looks like this

```shell
ls -l ~/tools/
total 861M
drwxr-xr-x  7 lie128 hpc-users 4.0K Sep 22 13:19 megalodon/
drwxr-xr-x 10 lie128 hpc-users 8.0K Sep 22 13:07 minimap2/
drwxr-xr-x  5 lie128 hpc-users 4.0K Sep 22 11:39 ont-guppy/
-rw-r--r--  1 lie128 hpc-users 857M Sep 22 10:53 ont-guppy_5.0.14_linux64.tar.gz
drwxr-xr-x  5 lie128 hpc-users 4.0K Sep 22 11:39 rerio/
```

## Troubleshooting corner ##

We're still not 100% sure WHY these errors occur, but google + elbow grease led to fixes that work (mysteriously). From brief discussion with CSIRO Scientific Computing, oddities could sometimes occur when mixing `conda` packages (usually compiled with `gcc`) and `pip install` stuff (which sometimes uses the system default compiler, `icc` in this case).

Problem 1: `mappy` (Python bindings for `minimap2`, automatically installed by `pip install megalodon`) complains about intel sse2 something something.

```py
Traceback (most recent call last):
  File "/home/lie128/.conda/envs/megalodon/bin/megalodon", line 8, in <module>
    sys.exit(_main())
  File "/home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/__main__.py", line 689, in _main
    from megalodon import megalodon
  File "/home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/megalodon.py", line 11, in <module>
    import mappy
ImportError: /home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/mappy.cpython-38-x86_64-linux-gnu.so: undefined symbol: __intel_sse2_strcpy
```

Fix: recompile `mappy` with `gcc`

```shell
pip remove mappy
module load gcc/11.1.0
export CC=gcc
cd ~/tools
git clone https://github.com/lh3/minimap2
cd minimap2
pip install .
```

Problem 2: `megalodon` complains about... svml something something.

```py
Traceback (most recent call last):
  File "/home/lie128/.conda/envs/megalodon/bin/megalodon", line 8, in <module>
    sys.exit(_main())
  File "/home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/__main__.py", line 689, in _main
    from megalodon import megalodon
  File "/home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/megalodon.py", line 15, in <module>
    from megalodon import (
  File "/home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/aggregate.py", line 9, in <module>
    from megalodon import (
  File "/home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/mods.py", line 15, in <module>
    from megalodon import (
ImportError: /home/lie128/.conda/envs/megalodon/lib/python3.8/site-packages/megalodon/decode.cpython-38-x86_64-linux-gnu.so: undefined symbol: __svml_log1pf4
```

Fix: also recompile `megalodon` with `gcc`

```shell
pip remove megalodon
module load gcc/11.1.0
export CC=gcc
cd ~/tools
git clone https://github.com/nanoporetech/megalodon
cd megalodon
pip install .
```
