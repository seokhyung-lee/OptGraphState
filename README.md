# OptGraphState

**Version 0.1.1**

*Graph-theoretical optimization of fusion-based graph state generation.*

**OptGraphState** is a python package that implements the graph-theoretical strategy to optimize the fusion-based
generation of any graph state. The package has the following features:

- Finding a resource-efficient method of generating a given graph state through type-II fusions from multiple basic
  resource states, which are three-qubit linear graph states.
- Computing the corresponding resource overhead, which is quantified by the average number of required basic resource
  states.
- Visualizing the original graph (of the graph state you want to generate), unraveled graphs, and fusion networks. An
  unraveled graph is a simplified graph where the corresponding graph state is equivalent to the desired graph state up
  to fusions and single-qubit Clifford operations. A fusion network is a graph that instructs the fusions between basic
  resource states required to generate the desired graph state.
- Various predefined sample graphs for input.

Tutorials: https://github.com/seokhyung-lee/OptGraphState/raw/main/tutorials.pdf

API reference: https://seokhyung-lee.github.io/OptGraphState