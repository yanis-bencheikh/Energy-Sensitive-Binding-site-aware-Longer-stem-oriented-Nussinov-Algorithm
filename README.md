# Enhanced Nussinov Algorithm for RNA Secondary Structure Prediction

An energy-sensitive, binding-site-aware, and longer-stem-oriented implementation of the Nussinov algorithm with minimal loop length thresholding for accurate RNA secondary structure prediction.

## 🧬 Overview

This project presents novel variants of the classical Nussinov-Jacobson algorithm that significantly improve RNA secondary structure prediction accuracy. Our best-performing model achieves **89% accuracy** - a **threefold improvement** over the classical Nussinov algorithm while maintaining the same O(n³) computational complexity.

## ✨ Key Features

### 🎯 Enhanced Algorithms
- **BS&E Model**: Binding-site aware and energy-sensitive prediction
- **BS&E+AltTB Model**: Includes alternative traceback for longer stems and loops
- **NJMLL Model**: Classical Nussinov with minimal loop length constraints

### 🔬 Biological Constraints Integration
- **Minimal Loop Length**: Prevents biologically unstable short loops (< 3 bases)
- **Energy-Sensitive Base Pairing**: Uses Clote and Backofen scoring scheme
- **Binding Site Awareness**: Excludes known binding sites from base pairing
- **Stem Length Prioritization**: Alternative traceback favors longer, more stable stems

### 📊 Performance Highlights
- **89% average accuracy** (BS&E+AltTB model)
- **Perfect prediction (100% accuracy)** for 6 out of 12 test RNA sequences
- **249% improvement** over classical Nussinov algorithm
- **O(n³) runtime complexity** maintained

## 🚀 Getting Started

### Prerequisites
```bash
python >= 3.7
numpy
pandas
matplotlib (for visualization)
```

### Installation
```bash
git clone https://github.com/yourusername/enhanced-nussinov
cd enhanced-nussinov
pip install -r requirements.txt
```

### Basic Usage

```python
# Import the enhanced algorithms
from enhanced_nussinov import binding_energy_and_alternate_tb_nussinov

# Your RNA sequence
rna_sequence = "GGGAAAUCC"

# Option 1: Best performing model (BS&E+AltTB)
predictor = binding_energy_and_alternate_tb_nussinov()
base_pairs = predictor.run(rna_sequence, ["AUC"], verbose=False, minimal_loop_length=5)

# Option 2: Optimize minimal loop length
def find_optimal_mll(rna, binding_sites, max_length=10):
    best_accuracy = 0
    optimal_mll = 0
    
    for mll in range(max_length):
        # Test with different lengths
        predictor = binding_energy_and_alternate_tb_nussinov()
        result = predictor.run(rna, binding_sites, False, mll)
        # Evaluate accuracy (requires observed structure for comparison)
        
    return optimal_mll

# Option 3: Visualize results
def visualize_structure(rna, base_pairs):
    dot_notation = dot_write(rna, base_pairs)
    print(f"RNA Sequence: {rna}")
    print(f"Structure:    {dot_notation}")
    return dot_notation
```

## 📋 Available Models

### 1. Classical Nussinov (NJ)
Basic implementation with base pair maximization.
```python
from enhanced_nussinov import basic_nussinov
predictor = basic_nussinov()
result = predictor.run(rna_sequence, verbose=False, minimal_loop_length=0)
```

### 2. Nussinov with Minimal Loop Length (NJMLL)
Adds minimal loop length constraint.
```python
from enhanced_nussinov import basic_nussinov
predictor = basic_nussinov()
result = predictor.run(rna_sequence, verbose=False, minimal_loop_length=5)
```

### 3. Binding Site & Energy Aware (BS&E)
Incorporates energy model and binding site awareness.
```python
from enhanced_nussinov import binding_and_energy_nussinov
predictor = binding_and_energy_nussinov()
result = predictor.run(rna_sequence, ["AUC", "GGG"], verbose=False, minimal_loop_length=5)
```

### 4. BS&E with Alternative Traceback (BS&E+AltTB)
**Recommended**: Best performing model with stem length prioritization.
```python
from enhanced_nussinov import binding_energy_and_alternate_tb_nussinov
predictor = binding_energy_and_alternate_tb_nussinov()
result = predictor.run(rna_sequence, ["AUC"], verbose=False, minimal_loop_length=5)
```

## 🧪 Benchmarking and Validation

### Test Dataset
Our algorithms were benchmarked on 12 RNA sequences including:
- Tar 16, Tar 16*
- R1inv, R2inv  
- DIS, CopA, CopT
- ATP Sensitive Ribozyme
- lncRNA54, RepZ
- RyhB, SodB, OxyS, fhlA

### Performance Comparison

| Model | Mean Accuracy | Improvement vs NJ |
|-------|---------------|-------------------|
| NJ (Classical) | 25.6% | - |
| NJMLL | 78.2% | +205% |
| BS&E | 82.8% | +223% |
| **BS&E+AltTB** | **89.4%** | **+249%** |

### Perfect Predictions (100% Accuracy)
The BS&E+AltTB model achieved perfect accuracy for:
- Tar 16, Tar 16*
- R1inv, R2inv
- CopT, lncRNA54

## 🔧 Advanced Usage

### Hyperparameter Optimization
```python
import math
import matplotlib.pyplot as plt

def minimal_loop_length_effect(rna, rna_name, observed_structure, binding_sites):
    """Optimize minimal loop length for a given RNA sequence"""
    accuracy_dict = {}
    max_accuracy = 0
    optimal_mll = 0
    
    # Test different loop lengths
    for mll in range(round(math.sqrt(len(rna)))):
        # Evaluate all models
        accuracies = evaluate_acc_all_models(rna, observed_structure, mll, binding_sites, False)
        accuracy_dict[f'{mll}'] = accuracies
        
        # Find best accuracy
        for acc in accuracies:
            if acc > max_accuracy:
                max_accuracy = acc
                optimal_mll = mll
    
    # Visualize results
    plt.figure(figsize=(10, 6))
    plt.title(f'Minimal Loop Length Optimization - RNA: {rna_name}')
    plt.axvline(x=optimal_mll, color='purple', ls='--', label=f'Optimal MLL: {optimal_mll}')
    plt.xlabel("Minimal loop length")
    plt.ylabel("Prediction accuracy")
    
    models = ['basic', 'basic_MLL', 'binding_energy_MLL', 'binding_energy_altTB_MLL']
    for i, model in enumerate(models):
        accuracies = [accuracy_dict[str(mll)][i] for mll in range(len(accuracy_dict))]
        plt.plot(range(len(accuracy_dict)), accuracies, '-', label=model)
    
    plt.legend()
    plt.show()
    
    return optimal_mll

# Example usage
rna = "AAGACACUCCGGGGGUGAUAGAAAGCAGCAAGGCGGUUCAAGCAUUCUUUCUAUCACCCCCAAAAGGAAAAUACCG"
observed = get_pairs_from_dot(".......((((((((((((((((((..((..((....))..))...)))))))))))))))....)))........")
binding_sites = ["AAGCAGCAAGGC", "UCAAGCAUUCUU"]

optimal_mll = minimal_loop_length_effect(rna, "RepZ", observed, binding_sites)
```

### Accuracy Evaluation
```python
def evaluate_acc(observed_structure, predicted_structure, verbose=False):
    """Evaluate base pair prediction accuracy"""
    obs_ord = sorted(observed_structure)
    pred_ord = sorted(predicted_structure)
    
    correct = 0
    for bp in pred_ord:
        if bp in obs_ord:
            correct += 1
    
    accuracy = correct / len(obs_ord) if len(obs_ord) > 0 else 0
    
    if verbose:
        print(f'Predicted pairs: {pred_ord}')
        print(f'Observed pairs: {obs_ord}')
        print(f'Accuracy: {accuracy:.3f}')
    
    return accuracy

def get_pairs_from_dot(dot_notation):
    """Convert dot-bracket notation to base pairs"""
    stack = []
    pairings = []
    
    for i, char in enumerate(dot_notation):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening_index = stack.pop()
                pairings.append((opening_index, i))
    
    return pairings
```

### Structure Visualization
```python
def dot_write(rna, fold, binding_sites_indices=None):
    """Generate dot-bracket notation from base pairs"""
    dot = ["."] * len(rna)
    
    # Mark binding sites with asterisks
    if binding_sites_indices:
        for i in binding_sites_indices:
            dot[i] = "*"
    
    # Mark base pairs
    for pair in fold:
        dot[min(pair)] = "("
        dot[max(pair)] = ")"
    
    return "".join(dot)

def compare_structures(rna, predicted_pairs, observed_pairs, binding_sites=None):
    """Compare predicted and observed structures"""
    print(f"RNA Sequence: {rna}")
    print(f"Predicted   : {dot_write(rna, predicted_pairs, binding_sites)}")
    print(f"Observed    : {dot_write(rna, observed_pairs, binding_sites)}")
    
    accuracy = evaluate_acc(observed_pairs, predicted_pairs, verbose=True)
    return accuracy
```

## 🧮 Energy Model

Our energy-sensitive models use the Clote and Backofen scoring scheme:

| Base Pair | Energy Value |
|-----------|--------------|
| G-C | -5 |
| A-U | -4 |
| G-U | -1 |

This accounts for the varying stability of different base pairings, with G-C pairs being most stable.

## 📊 Minimal Loop Length Optimization

The algorithms support automatic optimization of minimal loop length:

```python
def find_optimal_mll(rna_sequence, max_mll=None):
    """Find optimal minimal loop length"""
    if max_mll is None:
        max_mll = int(len(rna_sequence)**0.5)
    
    results = []
    for mll in range(0, max_mll):
        # Test with this MLL
        predictor = binding_energy_and_alternate_tb_nussinov()
        # Evaluate performance (requires reference structure)
        accuracy = evaluate_performance(rna_sequence, mll)
        results.append((mll, accuracy))
    
    optimal_mll = max(results, key=lambda x: x[1])[0]
    return optimal_mll, results
```

## 🔬 Biological Applications

### Long Non-coding RNAs (lncRNAs)
Our algorithms are particularly effective for:
- Gene expression regulation prediction
- Chromatin function analysis
- mRNA stability assessment
- Signaling pathway modeling

### RNA-RNA Interactions
- Antisense RNA binding prediction
- miRNA target prediction (with binding site constraints)
- Ribozyme structure analysis

## ⚡ Performance Characteristics

- **Time Complexity**: O(n³) where n is sequence length
- **Space Complexity**: O(n²) for dynamic programming table
- **Scalability**: Tested on sequences up to 100+ nucleotides
- **Memory Efficient**: Minimal overhead over classical Nussinov

## 📁 Project Structure

```
enhanced-nussinov/
├── copie_de_project_lncrna_comp561.py    # Main source code
├── README.md                              # This file
├── requirements.txt                       # Python dependencies
├── data/                                  # Test data
│   ├── test_sequences.txt
│   └── observed_structures.txt
├── examples/                              # Usage examples
│   ├── basic_usage.py
│   ├── optimization_example.py
│   └── visualization_demo.py
└── docs/                                  # Documentation
    ├── algorithm_details.md
    └── performance_analysis.md
```

## 🔄 Evaluation Pipeline

```python
def evaluate_acc_all_models(rna, observed_structure, minimal_loop_length, binding_sites, verbose=False):
    """Evaluate all models and compare their performance"""
    
    # Initialize models
    m1 = basic_nussinov()
    m2 = binding_and_energy_nussinov()
    m3 = binding_energy_and_alternate_tb_nussinov()
    
    # Get predictions
    pred0 = m1.run(rna, False, 0)  # Classical NJ
    pred1 = m1.run(rna, False, minimal_loop_length)  # NJMLL
    pred2 = m2.run(rna, binding_sites, False, minimal_loop_length)  # BS&E
    pred3 = m3.run(rna, binding_sites, False, minimal_loop_length)  # BS&E+AltTB
    
    # Calculate accuracies
    acc0 = evaluate_acc(observed_structure, pred0, verbose)
    acc1 = evaluate_acc(observed_structure, pred1, verbose)
    acc2 = evaluate_acc(observed_structure, pred2, verbose)
    acc3 = evaluate_acc(observed_structure, pred3, verbose)
    
    if verbose:
        print(f'\nClassical Nussinov: {acc0:.3f}')
        print(f'Nussinov MLL: {acc1:.3f}')
        print(f'BS&E: {acc2:.3f}')
        print(f'BS&E+AltTB: {acc3:.3f}')
    
    return [acc0, acc1, acc2, acc3]
```

## 🤝 Contributing

We welcome contributions with enthusiasm! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Areas for Enhancement
- Pseudoknot prediction integration
- Additional thermodynamic parameters
- GPU acceleration
- Ensemble prediction methods
- Graphical user interface
- Support for additional file formats (FASTA, etc.)

### How to Contribute
1. **Fork** the project
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a **Pull Request**

## 📖 Citation

If you use this work in your research, please cite:

```bibtex
@article{cao2023enhanced,
    title={An Energy Sensitive, Binding-site-aware and Longer-stem-oriented Nussinov Algorithm with Minimal Loop Length Thresholding},
    author={Cao, Sophie and Bencheikh, Yanis and Lara, Vitoria},
    journal={COMP561 Final Project},
    year={2023},
    publisher={McGill University},
    note={Enhanced Nussinov algorithm for RNA secondary structure prediction}
}
```

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Based on the classical Nussinov-Jacobson algorithm (1980)
- Energy scoring adapted from Clote and Backofen methodology
- Inspired by RNA folding research in molecular biology
- Special thanks to the COMP561 course at McGill University
- Thanks to the RNA research community for valuable feedback

## 🔬 Authors

- **Sophie Cao** - Implementation of BS&E and BS&E+AltTB models, report sections
- **Yanis Bencheikh** - Implementation contribution, MLL optimization framework, algorithmic pipeline
- **Vitoria Lara Soria** - Literature review, RNA sequence annotation, test dataset composition

## 📞 Support

For questions, issues, or collaboration opportunities:
- 📧 Email: [your-email@domain.com]
- 🐛 Issues: [GitHub Issues](https://github.com/yourusername/enhanced-nussinov/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/yourusername/enhanced-nussinov/discussions)
- 📚 Documentation: [Project Wiki](https://github.com/yourusername/enhanced-nussinov/wiki)

## 🎯 Roadmap

### Version 2.0 (Planned)
- [ ] Pseudoknot support
- [ ] Graphical user interface
- [ ] REST API for online predictions
- [ ] Additional file format support
- [ ] Performance optimizations

### Version 1.5 (In Progress)
- [x] Automatic hyperparameter optimization
- [x] Enhanced visualizations
- [ ] Interactive documentation
- [ ] Comprehensive unit tests

---

**Made with ❤️ for the RNA research community** 🧬✨

*"Understanding RNA structure is understanding life itself"*
