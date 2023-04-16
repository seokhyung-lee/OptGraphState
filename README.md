# OptGraphState

*Graph-theoretical optimization of fusion-based graph state generation.*

**OptGraphState** is a python package that implements the graph-theoretical
strategy to optimize the fusion-based generation of any graph state.
The package has the following features:

- Implementation of individual steps of the strategy (unraveling a given graph,
constructing a fusion network, and determining the fusion order) and
computation of the corresponding resource overhead (number of three-qubit
linear graph states required to generate the graph state).
- Iteration of these three steps (with a fixed iteration number or the adaptive
iteration method) for minimizing the resource overhead. Multiprocessing is
natively supported.
- Visualization of the original graph, unraveled graph, and fusion network.
- Provision of a clear instruction for building a desired graph state from
three-qubit linear graph states.
- Various predefined sample graphs for input.

Tutorials: https://github.com/seokhyung-lee/OptGraphState/raw/main/tutorials.pdf

API reference: https://seokhyung-lee.github.io/OptGraphState
