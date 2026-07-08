## 04_Classification_lncRNAs_by_Genomic_Location

### Purpose

Classify identified lncRNAs according to their genomic location relative to neighboring protein-coding genes (e.g., intergenic, antisense, intronic, sense-overlapping).

### Tool

[FEELnc (FlExible Extraction of LNcRNA)](https://github.com/tderrien/feelnc) v0.2.1

> Wucher, V. et al. (2017). *FEELnc: a tool for long non-coding RNA annotation and its application to the dog transcriptome.* Nucleic Acids Research. <https://doi.org/10.1093/nar/gkw1306>

### Platform

Local (terminal)

### Input / Output

|  |  |
|----|----|
| Input | lncRNA GTF, protein-coding (mRNA) reference GTF |
| Output | Genomic lncRNA classification (interaction/class table), Transcripts with isBest? == 1 were retained. |

------------------------------------------------------------------------

### Installation

FEELnc is installed via a dedicated `conda` environment (bioconda). On **Apple Silicon (M1/M2/M3) Macs**, some dependencies (e.g., `r-base=3.2.2`, older BioConda/BioBuilds packages) are only available as `x86_64` builds, so an Intel-emulated (Rosetta 2) shell and Miniconda installation are required first.

> For installation problems specific to your OS/architecture, please refer to the [GitHub repository](https://github.com/tderrien/feelnc).

#### Step 1 — Enable Rosetta 2 and open an x86_64 shell (macOS Apple Silicon only)

``` bash
# Install Rosetta 2 (one-time setup)
softwareupdate --install-rosetta

# Close and reopen the terminal, then launch an x86_64-emulated shell
arch -x86_64 zsh
```

#### Step 2 — Install Miniconda (x86_64 build)

``` bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh -u
```

#### Step 3 — Install Perl and required CPAN dependencies

``` bash
brew install perl
curl -L https://cpanmin.us | perl - App::cpanminus

# Locate the installed cpanm binary
find /usr/local ~/perl5 -type f -name cpanm 2>/dev/null

# Add Perl tool paths to the environment
export PATH="$HOME/perl5/bin:$PATH"
export PERL5LIB="$HOME/perl5/lib/perl5${PERL5LIB:+:$PERL5LIB}"
```

> Add the two `export` lines to your shell profile (`~/.zshrc` or `~/.bash_profile`) to persist them across sessions.

#### Step 4 — Configure conda channels

``` bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels BioBuilds
conda config --set channel_priority strict
```

> If dependency resolution fails under `strict` priority, fall back to:
>
> ``` bash
> conda config --set channel_priority flexible
> ```

#### Step 5 — Create and activate the FEELnc environment

``` bash
conda create -n feelnc_env \
  -c conda-forge -c BioBuilds \
  perl-bioperl=1.6.924 \
  perl-parallel-forkmanager \
  r-base=3.2.2 \
  bedtools \
  mamba

conda activate feelnc_env
```

#### Step 6 — Install FEELnc

``` bash
conda install -c bioconda feelnc --no-deps
```

#### Verify installation

``` bash
feelnc_classifier.pl --help
```

------------------------------------------------------------------------

### Usage

``` bash
FEELnc_classifier.pl -i <lncRNA>.gtf -a <mRNA_reference>.gtf > <output_dir>/lncRNA_classes.txt
```

**Example runs** (classifying different lncRNA sets against the same araport11 reference annotation):

``` bash
FEELnc_classifier.pl -i lncRNAs_annotated.gtf -a araport11.gtf > Results/lncRNA_annotated_classes.txt
FEELnc_classifier.pl -i lncRNAs_novel.gtf     -a araport11.gtf > Results/lncRNA_novel_classes.txt
```

### Options

**General**

| Flag              | Description                                        |
|-------------------|----------------------------------------------------|
| `-b, --biotype`   | Print the biotype of each transcript in the output |
| `-l, --log=file`  | Specify the name of the log file                   |
| `-v, --verbosity` | Level of verbosity                                 |
| `--help`          | Print help message                                 |
| `--man`           | Open the man page                                  |

**Mandatory arguments**

| Flag                    | Description                        |
|-------------------------|------------------------------------|
| `-i, --lncrna=file.gtf` | lncRNA GTF file                    |
| `-a, --mrna=file.gtf`   | Protein-coding annotation GTF file |

**Filtering arguments**

| Flag | Description |
|----|----|
| `-w, --window=10000` | Window size around the lncRNA for interaction/classification (default: 10,000 bp) |
| `-m, --maxwindow=100000` | Maximum window size during the expansion process (default: 100,000 bp) |

------------------------------------------------------------------------
