# MaximumBiclique
An implementation of the paper "Theoretically and Practically Efficient Maximum Biclique Search"

# Compile
Using cmake files to compile the source files. For instance,
```
cmake .
make
```

# Usage
To execute the code, you need to run the following executable files, which accept the following optional parameters:

- "-f": the running bipartite graph.

- "-l": The left size constraint.

- "-r": The right size constraint.

- "-a": the algorithm to run ("-a 1" for the algorithm MCES, "-a 2" for the algorithm iMCES). 

An running example:

```
.bin/run -f datas/pics.txt -l 5 -r 5 -a 2
```
