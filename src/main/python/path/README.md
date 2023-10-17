## Predicting Affinity Through Homology (PATH): Interpretable Binding Affinity Prediction with Persistent Homology

*This repository contains the optimized inference code only*

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

### Running the inference algorithm
`python inference.py`
