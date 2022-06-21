<p align="center"><img align="center" width="217" height="175" src="https://github.com/gwct/bonsai/blob/main/img/bonsai.png?raw=true"></p>

# Bonsai

## Tree pruning using concordance factors and branch lengths

# Authors

### Gregg Thomas and Jeffrey Good

## About

Bonsai offers informed tree pruning to maximize the phylogenetic concordance of the loci underlying a species tree. This can help facilitate analyses on large trees.

The program also offers a suite of alignment and concordance factor calculations.

## Installation

Simply download the program and run it. You may want to add the Referee folder to your $PATH variable for ease of use.
### The only dependency is Python 3+

## Usage

1. Label a species tree:

```bash
python bonssai.py -st [species tree file] --labeltree
```

2. Calculate concordance factors using input gene tree and exit:

```bash
python bonssai.py -st [species tree file] -gt [gene trees file] --cf
```

3. Calculate concordance factors using input gene tree and alignments, and calculate alignment statistics and exit:

```
python bonssai.py -st [species tree file] -gt [gene trees file] -d [directory with alignments in FASTA format] --stats --cf
```

3. Use input gene trees to calculate concordance factors and prune the species tree based on them. Do a maximum of 3 iterations of pruning:

```
python bonssai.py -st [species tree file] -gt [gene trees file] -i 3
```

3. Prune the input species tree using the labels already in the tree tree based on them. Do a maximum of 5 iterations of pruning:

```
python bonssai.py -st [species tree file] -gt [gene trees file] -i 5 --labels
```

## Options

| Option | Description | 
| ------ | ----------- |
| `-st` | A file or string containing a Newick formatted tree as the species tree (REQUIRED). |
| `-gt` | A file containing one or more Newick formatted trees, one per line, as gene trees (REQUIRED except with `--labeltree` or `--prune`). |
| `-o` | The desired output directory. This will be created for you if it doesn't exist. Default: `bonsai-[date]-[time]` |
| `-b` | The lower percentile of branch lengths to consider for paring at each iteration. Must be between 0 and 100. Default: 10 |
| `-g` | The lower threshold of gene concordance factor to consider for paring at each iteration. Must be between 0 and 100. Default: 25. |
| `-i` | The maximum number of paring iterations to perform before stopping. Default: 3. |
| `-m` | The maximum number of species (tips) that can be pruned while paring any given single branch. Default: 10. |
| `-p` | A file containing branches to prune from the input tree(s), one branch per line defined as internal node labels with `--labeltree` or 2 tips that descend from that branch (e.g. 'spec1 spec2'). Lines in file starting with '#' are ignored. |
| `-e` | A file containing branches to be exempt from pruning (note that they can still be PARED), one branch per line defined as internal node labels with `--labeltree` or 2 tips that descend from that branch (e.g. 'spec1 spec2'). Lines in file starting with '#' are ignored |
| `-d` | A directory of corresponding alignments to calculate site concordance factors on the species tree (optional). |
| `-scf` | The number of quartets to sample around each branch for sCF calculations. |
| `-n` | The number of processes that Bonsai should use. Default: 1. |
| `--cf` | Stop after calculating concordance factors. |
| `--stats` | When a directory of alignments is given with `-scf`, set this flag to write out some stats for the alignments. |
| `--labels` | Set this option to disable concordance factor calculations and use the labels on the input tree for paring. Note that if these labels are NOT concordance factors, you may also need to adjust `-g`. |
| `--labeltree` | Just read the species tree from the input, label the internal nodes, print, and exit. |
| `--overwrite` | Set this to overwrite existing files. |
| `--appendlog` | Set this to keep the old log file even if `--overwrite` is specified. New log information will instead be appended to the previous log file. |
| `--version` | Simply print the version and exit. Can also be called as `-version`, `-v`, or `--v`. |
| `--info` |  Print some meta information about the program and exit. No other options required. |
| `--quiet` | Set this flag to prevent degenotate from reporting detailed information about each step. |



