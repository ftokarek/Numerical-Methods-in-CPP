# Cholesky Decomposition Implementation

## Overview

This project implements **Cholesky decomposition** in C++ using modern, professional techniques. Cholesky decomposition is a specialized factorization method for symmetric positive definite (SPD) matrices, decomposing A into A = L·L^T where L is a lower triangular matrix with positive diagonal elements. 
This method is approximately **twice as fast** as Gaussian elimination for SPD matrices.

## Mathematical Background

### Cholesky Decomposition Theory

For a symmetric positive definite matrix A, Cholesky decomposition finds a unique lower triangular matrix L such that:

```
A = L · L^T
```

Where:
- **L** is lower triangular with positive diagonal elements
- **L^T** is the transpose of L (upper triangular)
- The decomposition exists **if and only if** A is symmetric positive definite

### Algorithm Steps

1. **Decomposition**: Compute L such that A = L·L^T
2. **Forward Substitution**: Solve L·y = b for y
3. **Backward Substitution**: Solve L^T·x = y for x

### Mathematical Formulation

The elements of L are computed as:

```
for j = 0 to n-1:
    for i = j to n-1:
        if i == j:
            L[i][j] = √(A[i][i] - Σ(k=0 to j-1) L[i][k]²)
        else:
            L[i][j] = (A[i][j] - Σ(k=0 to j-1) L[i][k]·L[j][k]) / L[j][j]
```

## Features

- **Optimal Performance**: ~2x faster than Gaussian elimination for SPD matrices
- **Numerical Stability**: Specialized algorithm for symmetric positive definite matrices
- **Comprehensive Validation**: 
  - Symmetry checking
  - Positive definiteness verification
  - Input validation and error handling
- **Professional C++ Design**:
  - Modern C++17 features
  - Exception safety and const correctness
  - Comprehensive documentation
  - RAII and move semantics
- **Solution Verification**: Residual computation for accuracy checking
- **Matrix Properties**: Determinant calculation and condition assessment
- **Memory Efficient**: Optimized storage for triangular matrices

## Files Structure

```
04_cholesky_decomposition/
├── README.md                    # This file
├── main.cpp                     # Main program with file input
├── cholesky_solver.hpp          # Header with CholeskySolver class
├── data.txt                     # Sample input data file
└── main                         # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::CholeskySolver`
```cpp
class CholeskySolver 
{
    public:
        explicit CholeskySolver(std::vector<std::vector<double>> matrix);
        [[nodiscard]] std::vector<double> solve(const std::vector<double>& b);
        [[nodiscard]] const std::vector<std::vector<double>>& get_L_matrix() const;
        [[nodiscard]] double verify_solution(const std::vector<double>& b) const;
        [[nodiscard]] double get_determinant() const;
        [[nodiscard]] bool is_positive_definite() const;
        [[nodiscard]] std::size_t size() const noexcept;
};
```

## Usage Example

```cpp
#include "cholesky_solver.hpp"
using namespace numerical_methods;

// Create symmetric positive definite matrix
std::vector<std::vector<double>> matrix = {
    {4,  12, -16},
    {12, 37, -43},
    {-16, -43, 98}
};

std::vector<double> rhs = {1, 2, 3};

// Create solver and solve
CholeskySolver solver(std::move(matrix));
auto solution = solver.solve(rhs);

// Verify solution
double residual = solver.verify_solution(rhs);
```

## Input Data Format

The program reads data from `data.txt` with the following format:

```
n
a₁₁ a₁₂ ... a₁ₙ
a₂₁ a₂₂ ... a₂ₙ
...
aₙ₁ aₙ₂ ... aₙₙ
b₁ b₂ ... bₙ
```

### Example `data.txt`:
```
3
4 12 -16
12 37 -43
-16 -43 98
1 2 3
```

## Compilation

```bash
# Compile the program
g++ -std=c++17 -O3 -Wall -Wextra main.cpp -o main
```

```bash
# Run the program
./main
```

or

```bash
# Compile the program
clang++ -std=c++17 main.cpp -o main
```

```bash
# Run the program
./main
```

## Sample Output

```
Cholesky Decomposition Solver
============================================================
Successfully loaded 3x3 system from data.txt

Coefficient Matrix A:
--------------------------------------------------
[      4.000     12.000    -16.000 ]
[     12.000     37.000    -43.000 ]
[    -16.000    -43.000     98.000 ]
--------------------------------------------------

============================================================
CHOLESKY DECOMPOSITION RESULTS
============================================================

System size: 3x3

Lower triangular matrix L:
----------------------------------------
[    2.000    0.000    0.000 ]
[    6.000    1.000    0.000 ]
[   -8.000    5.000    3.000 ]

Solution:
------------------------------
x 1 =       28.583333
x 2 =       -7.666667
x 3 =        1.333333

Solution Verification:
------------------------------
Maximum residual: 2.664535e-14
✓ Solution is ACCURATE

Matrix Properties:
------------------------------
Determinant: 36.000000

============================================================
```

## Algorithm Complexity

- **Time Complexity**: O(n³/3) ≈ **half** of Gaussian elimination O(n³)
- **Space Complexity**: O(n²) for matrix storage
- **Numerical Stability**: Excellent for well-conditioned SPD matrices

## Requirements for Cholesky Decomposition

A matrix A must satisfy **both** conditions:

1. **Symmetric**: A = A^T (A[i][j] = A[j][i] for all i,j)
2. **Positive Definite**: All eigenvalues > 0 (equivalently, all leading principal minors > 0)

### Common SPD Matrix Examples:
- **Covariance matrices**
- **Gram matrices** (A^T·A where A has full rank)
- **Kernel matrices** in machine learning
- **Finite element stiffness matrices**

## Error Handling

The implementation provides comprehensive error handling:

- **Matrix Validation**: Symmetry and positive definiteness checking
- **Dimension Validation**: Square matrix and compatible RHS vector
- **Numerical Issues**: Detection of non-positive definite matrices
- **Input Validation**: File format checking, non-finite value detection

## Performance Advantages

Compared to Gaussian elimination:

| Feature | Cholesky | Gaussian |
|---------|----------|----------|
| **Operations** | n³/3 | 2n³/3 |
| **Speed** | ~2x faster | Baseline |
| **Memory** | n²/2 | n² |
| **Stability** | Excellent for SPD | Good with pivoting |
| **Applicability** | SPD matrices only | General matrices |

## Mathematical Properties

- **Uniqueness**: For SPD matrices, L is unique
- **Stability**: No pivoting needed for SPD matrices
- **Determinant**: det(A) = det(L)² = (∏ L[i][i])²
- **Condition**: Well-suited for well-conditioned SPD systems

## Applications

### Scientific Computing
- **Finite Element Methods**: Stiffness matrix solving
- **Optimization**: Quadratic programming, Newton's method
- **Statistics**: Multivariate normal distributions, regression

### Machine Learning
- **Gaussian Processes**: Kernel matrix inversion
- **Principal Component Analysis**: Covariance matrix decomposition
- **Neural Networks**: Second-order optimization methods

### Engineering
- **Structural Analysis**: Linear elastic systems
- **Signal Processing**: Autoregressive models
- **Control Theory**: Lyapunov equations

## Comparison with Other Methods

### vs. LU Decomposition
- **Faster**: Half the operations for SPD matrices
- **More Stable**: No pivoting required
- **Specialized**: Only works for SPD matrices

### vs. QR Decomposition
- **Faster**: O(n³/3) vs O(4n³/3)
- **Less General**: QR works for any matrix
- **Better Conditioned**: QR better for ill-conditioned systems

### vs. Gaussian Elimination
- **Efficiency**: 2x faster for SPD matrices
- **Stability**: Better numerical properties
- **Restrictions**: Requires SPD matrices

## Limitations

- **Matrix Requirements**: Only works for symmetric positive definite matrices
- **Numerical Precision**: Can fail for ill-conditioned matrices
- **Memory Usage**: Requires storing full L matrix
- **Input Validation**: Extensive checking needed for matrix properties

## Advanced Features

1. **Automatic Validation**: Checks symmetry and positive definiteness
2. **Solution Verification**: Computes residuals for accuracy assessment
3. **Determinant Calculation**: Efficient computation using L matrix
4. **Error Recovery**: Graceful handling of invalid inputs
5. **Memory Optimization**: Efficient storage patterns

## Future Enhancements

- **LDLT Decomposition**: For symmetric indefinite matrices
- **Incomplete Cholesky**: For sparse matrices
- **Parallel Cholesky**: Multi-threaded implementation
- **Iterative Refinement**: Improved accuracy for ill-conditioned systems
- **Band Cholesky**: Specialized for banded matrices

## Testing and Validation

The implementation includes comprehensive testing:
- **Known Solutions**: Verification against analytical solutions
- **Property Checking**: Symmetry and positive definiteness validation
- **Numerical Accuracy**: Residual analysis
- **Edge Cases**: Handling of near-singular matrices
