# COMBO: COmbing Multi-Bio Omics

[Ilaria Cosentini, Vincenza Barresi, Daniele Filippo Condorelli, Alfredo Ferro, Alfredo
Pulvirenti, and Salvatore Alaimo. Combo: A computational framework to
analyze rna-seq and methylation data through heterogeneous multi-layer networks.
In Hocine Cheri, Rosario Nunzio Mantegna, Luis M. Rocha, Chantal Cheri, and
Salvatore Micciche, editors, Complex Networks and Their Applications XI, pages
251{264, Cham, 2023. Springer International Publishing. ISBN 978-3-031-21127-0.](https://doi.org/10.1007/978-3-031-21127-0_21)

COMBO is a novel pipeline for multi-layer network inference and analysis to identify relevant pathways in the studied case. Taking advantage of the Boolean implication method, both transcriptomic and epigenomic data were analyzed through StepMiner[[1],[1]] and BooleanNet systems [[2],[2]]. The goal was to identify the implication between the different transcripts and methylated CpGs. The obtained results were used to generate heterogeneous multi-layer graphs. Subsequently, Neo4J was exploited to query the multi-layer network with properly defined Cypher queries.

#USAGE

```bash
bash COMBO.bash -h

Usage: COMBO.bash [options]

## Mandatory augument:
-E               expression matrix
-M               methylation matrix
-i               input file with sample information
-c               file with comparison between samples
-h               Print this Help
-o               output folder path

## Optional arguments:
# StepMiner and BooleanNet
-d <int>         delta threshold StepMiner (default=0.5)
-s <int>         statistic of an implication to be considered significant in BooleanNet (default=6.0)
-P <int>         maximum p-value P of an implication to be considered significant (default=0.01)

# Multilayer network creation
-w <string>      metapathway, can be KEGG or KEGG_Reactome (default=KEGG)
-f <int>         logFC threshold for expression data (default=|0.8|)
-x <int>         logFC threshold for methylation data (default=0)
-j <int>         adjusted pvalue threshold for methylation and expression data (default=0.05)
-n <string>      annotation name (mutation, protein expression)
-a               annotation table

# Number of CPU threads
-T <int>         threads number
```

# References
[[1]] Sahoo, D.,Dill, D.L., Tibshirani, R., Plevritis, S.K.: Extracting binary signals from microarray time-course data. Nucl. Acids Res., 3705–12 (2007)

[[2]] Sahoo, D., Dill, D.L., Gentles, A.J., Tibshirani, R., Plevritis, S.K.: Boolean implication networks derived from large scale, whole genome microarray datasets. Genome Biol. 9, R157 (2008)

[1]: https://academic.oup.com/nar/article/35/11/3705/2402546?login=false "Sahoo, D.,Dill, D.L., Tibshirani, R., Plevritis, S.K.: Extracting binary signals from microarray time-course data. Nucl. Acids Res., 3705–12 (2007)"

[2]: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-10-r157 "Sahoo, D., Dill, D.L., Gentles, A.J., Tibshirani, R., Plevritis, S.K.: Boolean implication networks derived from large scale, whole genome microarray datasets. Genome Biol. 9, R157 (2008)"
