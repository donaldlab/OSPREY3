## Predicting Affinity Through Homology (PATH): Interpretable Binding Affinity Prediction with Persistent Homology


Implementation of PATH+ and PATH- described in [*Predicting Affinity Through Homology (PATH): Interpretable Binding Affinity Prediction with Persistent Homology*](https://www.biorxiv.org/content/10.1101/2023.11.16.567384v3).

This repository contains the optimized inference code only. The training process and scripts for PATH are detailed [here](https://github.com/longyuxi/PATH-training). An open-source implementation of TNet-BP from TopologyNet (Cang and Wei, 2017) can be found [here](https://github.com/longyuxi/TopologyNet-2017).

### Prerequisites

**Required packages**: Python 3 (3.10 ≥ version ≥ 3.7), Numpy, Pandas 1.3.0, giotto-tda, Biopython, BioPandas, GUDHI

**Note:**
- Install giotto-tda using the method from [this issue](https://github.com/giotto-ai/giotto-tda/issues/671) if you want to use Python > 3.10.

**Installing these packages with conda:**

```
conda create -n PATH -c conda-forge python=3.10 scikit-learn=1.3.0 numpy pandas biopandas gudhi
conda activate PATH
pip install giotto-tda biopython pyfunctional
```

**Installing these packages without conda (make sure that you have 3.10 ≥ Python ≥ 3.7):**

```
pip install scikit-learn==1.3.0 numpy pandas biopandas gudhi giotto-tda biopython pyfunctional
```

### Running PATH+
`python path_plus/inference.py`

### Running PATH-
`python path_minus/pathminus_inference.py`

### Cite

```
@article{long2023predicting,
  title={Predicting Affinity Through Homology (PATH): Interpretable Binding Affinity Prediction with Persistent Homology},
  author={Long, Yuxi and Donald, Bruce},
  journal={bioRxiv},
  pages={2023--11},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```