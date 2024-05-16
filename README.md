# OptGraphState

**Version 0.3.1**

*Graph-theoretical optimization of fusion-based graph state generation.*

**OptGraphState** is a python package that implements the graph-theoretical strategy to optimize the fusion-based
generation of any graph state, which is proposed
in [Lee & Jeong, Quantum 7, 1212 (2023)](https://doi.org/10.22331/q-2023-12-20-1212).

The package has the following features:

- Finding a resource-efficient method of generating a given graph state through type-II fusions from multiple basic
  resource states, which are three-qubit linear graph states.
- Computing the corresponding resource overhead, which is quantified by the average number of required basic resource
  states or fusion attempts.
- Computing the success probability of graph state generation when the number of provided basic resource states is
  limited.
- Visualizing the original graph (of the graph state you want to generate), unraveled graphs, and fusion networks. An
  unraveled graph is a simplified graph where the corresponding graph state is equivalent to the desired graph state up
  to fusions and single-qubit Clifford operations. A fusion network is a graph that instructs the fusions between basic
  resource states required to generate the desired graph state.
- Various predefined sample graphs for input.

## Installation

`pip install optgraphstate`

## Prerequisites
The specified versions of these requirements are the ones with which we have tested OptGraphState. It is highly probable that newer versions work equally well.

(Exceptionally, we have checked that the latest version of `python-igraph` raises an error when visualizing graphs due to some internal changes of the library related to plotting graphs.)

- `python == 3.9`
- `numpy == 1.24.2`
- `python-igraph == 0.10.4`
- `networkx == 3.1`
- `matplotlib == 3.7.1`
- `parmap == 1.6.0` (Optional, for multiprocessing)
- `tqdm == 4.66.1` (Optiona, for progress bar in multiprocessing)

## Manuals

Tutorials: https://github.com/seokhyung-lee/OptGraphState/raw/main/tutorials.pdf

API reference: https://seokhyung-lee.github.io/OptGraphState

## License

OptGraphState is distributed under the MIT license. Please see the LICENSE file for more details.

## Citation

If you want to cite OptGraphState in an academic work, please cite the arXiv preprint:

```
@misc{lee2023graph,
      title={Graph-theoretical optimization of fusion-based graph state generation}, 
      author={Seok-Hyung Lee and Hyunseok Jeong},
      year={2023},
      eprint={2304.11988},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2304.11988}
}
```
