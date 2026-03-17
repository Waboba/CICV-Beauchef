# GRAIP: GRAphlet-based Incremental generator for Probabilistic networks

This is the official Python implementation of GRAIP, as described in  

Mornie B., Colle D., Audenaert P., Pickavet M. (2025) **Generating random graphs with prescribed graphlet frequency bounds derived from probabilistic networks.** PLOS ONE 20(8): e0328639. https://doi.org/10.1371/journal.pone.0328639

### Requirements  
To generate graphs with GRAIP, you need the following dependencies: (listed versions were used during development and testing)  
- Python: 3.12  
- numpy: 2.1.2
- matplotlib: 3.9.2  
- networkx: 3.3
- PyYAML: 6.0.2
  
Additionally, if you want to run `evaluate.py` and `test.py`, you also need:
- scikit-learn: 1.5.2

### Installation  
**1. Clone the repository**  
```bash
git clone https://github.com/bmornie/GRAIP.git
cd GRAIP
```
**2. Install required dependencies**  
```bash
pip install -r requirements.txt
```

### Test run
The `test.py` script provides an example of what can be done using our code. It includes: sampling, graph generation and evaluation of the generated graphs.
```bash
python test.py
```

### Basic usage
To generate graphs with GRAIP (or other models, see below) specify the parameters in `config.yaml` and run
```bash
python main.py
```
This only generates the graphs, but does not evaluate them. See `test.py` for an example on how to use the functions in `evaluate.py`, or read the function documentation.

### Code description
- `main.py` is the main script file. Parameters are specified in `config.yaml`.
- `evaluate.py` contains the functions used to perform the graph evaluation, as reported in the paper.
- `test.py` is a test script that shows how to use the GraphGenerator class and evaluation functions.
- `generator/generator.py` contains the GraphGenerator class. It is best to use this class instead of calling the functions in `graphlet_counts.py`, `sampling.py` and `models.py` directly. But running `main.py` should generally be sufficient.
  
In `models.py`, we also provide implementations of SwapCon and the (dual) BA model used for comparison in the paper. [GraphGen](https://github.com/idea-iitd/graphgen) is available on GitHub.

Finally, `data/` contains the PPI networks on which we evaluated our models. See `data/README.md` for more information on these datasets.
