# OptGraphState

**Version 0.1.2**

*Graph-theoretical optimization of fusion-based graph state generation.*

**OptGraphState** is a python package that implements the graph-theoretical strategy to optimize the fusion-based
generation of any graph state, which is proposed
in [Lee & Jeong, arXiv:2304.11988 [quant-ph] (2023)](https://arxiv.org/abs/2304.11988).

The package has the following features:

- Finding a resource-efficient method of generating a given graph state through type-II fusions from multiple basic
  resource states, which are three-qubit linear graph states.
- Computing the corresponding resource overhead, which is quantified by the average number of required basic resource
  states or fusion attempts.
- Visualizing the original graph (of the graph state you want to generate), unraveled graphs, and fusion networks. An
  unraveled graph is a simplified graph where the corresponding graph state is equivalent to the desired graph state up
  to fusions and single-qubit Clifford operations. A fusion network is a graph that instructs the fusions between basic
  resource states required to generate the desired graph state.
- Various predefined sample graphs for input.

## Installation

`pip install optgraphstate`

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