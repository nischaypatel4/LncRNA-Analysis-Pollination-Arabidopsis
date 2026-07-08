## 02_Coding_Potential_Prediction

### Purpose

Classify filtered transcripts as **coding** or **non-coding** based on sequence-intrinsic features, to identify candidate lncRNAs prior to downstream classification. Three complementary tools are used and their results are compared/intersected: **CPC2**, **CPAT**, and **LncFinder**. Transcripts retained as non-coding are further screened against the **Pfam-A** domain database to filter out any sequences with residual protein-coding signatures.

------------------------------------------------------------------------

## 1. CPC2

### Tool

[CPC2 (standalone)](https://github.com/gao-lab/CPC2_standalone) v2.0 Beta

> Kang, Y.-J., Yang, D.-C., Kong, L., Hou, M., Meng, Y.-Q., Wei, L., & Gao, G. (2017). *CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features.* Nucleic Acids Research, 45(W1), W12–W16. <https://doi.org/10.1093/nar/gkx428>

### Platform

Local (terminal)

### Input / Output

|        |                                                 |
|--------|-------------------------------------------------|
| Input  | Filtered FASTA (candidate transcript sequences) |
| Output | Coding/non-coding classification (text table)   |

### Installation

**Pre-requisite:** [Biopython](http://biopython.org/wiki/Download)

`CPC2-beta.tar.gz` can be downloaded from the [CPC2 download page](https://cpc2.gao-lab.org/download.php) and is also provided in this repository.

``` bash
# Unpack the tarball
gzip -dc CPC2-beta.tar.gz | tar xf -

# Build the bundled libsvm dependency
cd CPC2-beta
export CPC_HOME="$PWD"
cd libs/libsvm
gzip -dc libsvm-3.18.tar.gz | tar xf -
cd libsvm-3.18
make clean && make
```

> For system-specific installation problems, please refer to the [CPC2 download page](https://cpc2.gao-lab.org/download.php) or the [GitHub repository](https://github.com/gao-lab/CPC2_standalone).

### Usage

``` bash
cd $CPC_HOME
python ./bin/CPC2.py -i <input>.fa -o <output>.txt
```

------------------------------------------------------------------------

## 2. CPAT

### Tool

[CPAT](https://github.com/liguowang/cpat) v3.0.5

> Wang, L., Park, H. J., Dasari, S., Wang, S., Kocher, J.-P., & Li, W. (2013). *CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model.* Nucleic Acids Research, 41(6), e74. <https://doi.org/10.1093/nar/gkt006>

### Platform

Local (terminal)

### Input / Output

|        |                                                 |
|--------|-------------------------------------------------|
| Input  | Filtered FASTA (candidate transcript sequences) |
| Output | Coding probability scores                       |

### Installation

See the official [CPAT documentation](https://cpat.readthedocs.io/en/latest/) for installation specifications and requirements.

> For system-specific installation problems, please refer to the [GitHub repository](https://github.com/liguowang/cpat) or the associated [research paper](https://doi.org/10.1093/nar/gkt006).

### Usage

``` bash
cpat.py \
  -g <input>.fa \
  -d Plant.logit.RData \
  -x Plant_Hexamer.tsv \
  -o <output_prefix>
```

A coding probability **cutoff of 0.46** is applied to classify transcripts as coding vs. non-coding.

### Pre-built model / data provenance

The following plant-specific pre-built files are used and are provided in this repository

- `Plant_Hexamer.tsv`
- `Plant_model.rda`
- `Plant.logit.RData`
- Coding-probability cutoff: **0.46**

These were obtained from the **Plant-LncRNA-pipeline**:

> Tian, X.-C., Chen, Z.-Y., Nie, S., Shi, T.-L., Yan, X.-M., Bao, Y.-T., Li, Z.-C., Ma, H.-Y., Jia, K.-H., Zhao, W., & Mao, J.-F. (2024). *Plant-LncPipe: a computational pipeline providing significant improvement in plant lncRNA identification.* Horticulture Research, 11(4), uhae041. <https://doi.org/10.1093/hr/uhae041>
>
> Repository: <https://github.com/xuechantian/Plant-LncRNA-pipline>

------------------------------------------------------------------------

## 3. LncFinder

### Tool

[LncFinder](https://github.com/HAN-Siyu/LncFinder/) v1.1.6

> Han, S., Liang, Y., Ma, Q., Xu, Y., Zhang, Y., Du, W., Wang, C., & Li, Y. (2019). *LncFinder: an integrated platform for long non-coding RNA identification utilizing sequence intrinsic composition, structural information and physicochemical property.* Briefings in Bioinformatics, 20(6), 2009–2027. <https://doi.org/10.1093/bib/bby065>

### Platform

R / RStudio

### Input / Output

|        |                                                 |
|--------|-------------------------------------------------|
| Input  | Filtered FASTA (candidate transcript sequences) |
| Output | Coding/non-coding classification                |

### Installation

``` r
install.packages("LncFinder")
install.packages("seqinr")
```

> For system-specific installation problems, please refer to the [GitHub repository](https://github.com/HAN-Siyu/LncFinder/) or the associated [research paper](https://doi.org/10.1093/bib/bby065).

### Usage

``` r
library(LncFinder)
library(seqinr)

# Load reference training sequences (used to build k-mer frequency features)
mRNA   <- seqinr::read.fasta(file = "<path_to_data>/training_mRNA.fasta")
lncRNA <- seqinr::read.fasta(file = "<path_to_data>/training_lncRNA.fasta")

frequencies <- make_frequencies(
  cds.seq       = mRNA,
  lncRNA.seq    = lncRNA,
  SS.features   = FALSE,
  cds.format    = "DNA",
  lnc.format    = "DNA",
  check.cds     = TRUE,
  ignore.illegal = TRUE
)

# Load the pre-built plant SVM model
plant_model <- readRDS("<path_to_data>/Plant_model.rda")

# Load query sequences (candidate transcripts to classify)
query_seqs <- seqinr::read.fasta(file = "<path_to_input>/candidate_transcripts.fasta")

# Run classification
results <- LncFinder::lnc_finder(
  query_seqs,
  SS.features      = FALSE,
  format           = "DNA",
  frequencies.file = frequencies,
  svm.model        = plant_model,
  parallel.cores   = -1
)

# Write output
write.table(
  results,
  file      = "<path_to_output>/LncFinder_results.tsv",
  sep       = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote     = FALSE
)
```

A classification **cutoff of 0.5** is applied.

### Pre-built model / data provenance

The following plant-specific files are used and are provided in this repository

- `training_mRNA.fasta`
- `training_lncRNA.fasta`
- `Plant_model.rda`
- Classification cutoff: **0.5**

------------------------------------------------------------------------

## 4. Pfam (pfam_scan)

### Tool

[pfam_scan](https://github.com/aziele/pfam_scan) v1.6

> Refer to the [tool's repository](https://github.com/aziele/pfam_scan) for details. Domain annotations are based on the [Pfam-A](https://www.ebi.ac.uk/interpro/entry/pfam/) HMM database.

### Platform

Local (terminal)

### Database

Pfam-A HMM library (downloaded from the Pfam FTP site)

### Input / Output

|  |  |
|----|----|
| Input | Noncoding transcript candidates retained after CPC2, CPAT, and LncFinder filtering (translated to protein/ORF sequences in FASTA format) |
| Output | Pfam domain hits |

### Requirements

- Python \>= 3.8
- [HMMER](http://hmmer.wustl.edu/) \>= 3.3 (verify with `hmmscan -h`)

> For system-specific installation problems, please refer to the [GitHub repository](https://github.com/aziele/pfam_scan).

### Prepare the Pfam HMM library

``` bash
# 1. Download Pfam-A HMMs and their metadata
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# 2. Unpack into a dedicated database directory
mkdir pfamdb
gunzip -c Pfam-A.hmm.dat.gz > pfamdb/Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > pfamdb/Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz

# 3. Build HMMER binary index files
hmmpress pfamdb/Pfam-A.hmm
```

### Installation

`pfam_scan.py` is a single script and requires no formal installation — clone and run directly:

``` bash
git clone https://github.com/aziele/pfam_scan
cd pfam_scan
./pfam_scan.py --help
```

### Usage

``` bash
./pfam_scan.py <input>.fasta pfamdb/ -evalue 1e-3
```

Domain hits with an **e-value ≥ 1e−3 are discarded**, keeping only significant matches (`-evalue` sets the HMMER `--domE` domain-level e-value threshold).

------------------------------------------------------------------------
