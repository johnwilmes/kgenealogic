# kgenealogic

`kgenealogic` is a tool for genetic genealogy. It clusters relatives into a hypothesized family
tree based on shared and unshared DNA.

## Features

- Supports GEDmatch pairwise segment and triangulation files
- Uses both presence and inferred absence of triangulation
- Infers crossover probabilities on derived segments
- Allows, but doesn't require, specifying known relatives at each level of the family tree

## Usage

```
kgenealogic init [-p <project-file>]
kgenealogic add [-p <project-file>] ... # files to import
kgenealogic build [-p <project-file>] ... # build crossover probability model and find negative triangulations
kgenealogic cluster [-p <project-file>] [-o <out-file>] <config-file>
```
