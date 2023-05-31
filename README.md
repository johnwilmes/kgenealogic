# kgenealogic

`kgenealogic` is a tool for genetic genealogy. It clusters relatives into a hypothesized family
tree based on shared and unshared DNA.

## Features

- Supports GEDmatch pairwise segment and triangulation files
- Uses both presence and absence of triangulation
- Allows, but doesn't require, specifying known relatives at each node of the family tree
- Produces recommendations for additional GEDmatch exports that are most likely to improve accuracy
  of hypothesized tree

## Usage

```
kgenealogic init [-p <project-file>]
kgenealogic add [-p <project-file>] ... # files to import
kgenealogic build [-p <project-file>] ... # kits with complete triangulations
kgenealogic cluster [-p <project-file>] [-t <tree-file>] [-d <depth>] [-o <out-file>]
kgenealogic recommend
```
