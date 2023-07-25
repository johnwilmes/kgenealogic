# kgenealogic

`kgenealogic` is a tool for genetic genealogy. It clusters relatives into a hypothesized family
tree based on shared and unshared DNA.

## Features

- Supports GEDmatch pairwise segment and triangulation files
- Uses both presence and inferred absence of triangulation
- Allows, but doesn't require, specifying known relatives at each node of the family tree

## Usage

```
kgenealogic init [-p <project-file>]
kgenealogic add [-p <project-file>] ... # files to import
kgenealogic build [-p <project-file>] ... # kits with complete triangulations
kgenealogic cluster [-p <project-file>] [-d <depth>] [-o <out-file>] <tree-file>
```
