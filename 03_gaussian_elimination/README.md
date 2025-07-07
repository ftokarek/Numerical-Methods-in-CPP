# Gaussian Elimination Implementation

## Overview

This project implements **Gaussian elimination with partial pivoting** in C++ using modern, professional techniques. 
Gaussian elimination is a fundamental algorithm for solving systems of linear equations Ax = b by transforming the augmented matrix into row echelon form through forward elimination, followed by back substitution.

## Mathematical Background

Gaussian elimination transforms a system of linear equations:
```
a₁₁x₁ + a₁₂x₂ + ... + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + ... + a₂ₙxₙ = b₂
...
aₙ₁x₁ + aₙ₂x₂ + ... + aₙₙxₙ = bₙ
```

Into upper triangular form through **forward elimination**:
```
a'₁₁x₁ + a'₁₂x₂ + ... + a'₁ₙxₙ = b'₁
       a'₂₂x₂ + ... + a'₂ₙxₙ = b'₂
              ...
                    a'ₙₙxₙ = b'ₙ
```

Then solves for unknowns using **back substitution**.

### Partial Pivoting

To improve numerical stability, the algorithm uses **partial pivoting**:
- At each step, find the row with the largest absolute value in the current column
- Swap this row with the current pivot row
- This minimizes round-off errors and prevents division by small numbers

## Features

- **Numerical Stability**: Partial pivoting for robust computation
- **Comprehensive Error Handling**: Singular matrix detection and validation
- **Performance Optimized**: Efficient memory usage and algorithmic improvements
- **Solution Verification**: Residual computation for accuracy checking
- **Condition Number Estimation**: Numerical stability assessment
- **Determinant Calculation**: Additional matrix properties

## Files Structure

```
03_gaussian_elimination/
├── README.md                    # This file
├── main.cpp                     # Main program with file input
├── gaussian_solver.hpp          # Header with GaussianSolver class
├── data.txt                     # Sample input data file
└── main                         # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::GaussianSolver`
```cpp
class GaussianSolver 
{
    public:
        explicit GaussianSolver(std::vector<std::vector<double>> augmented_matrix);
        [[nodiscard]] std::vector<double> solve();
        [[nodiscard]] double verify_solution(const std::vector<std::vector<double>>& original_matrix,
                                            const std::vector<double>& rhs) const;
        [[nodiscard]] double get_condition_number() const;
        [[nodiscard]] double get_determinant() const;
        [[nodiscard]] bool is_well_conditioned() const;
        [[nodiscard]] std::size_t size() const noexcept;
};
```

## Usage Example

```cpp
#include "gaussian_solver.hpp"
using namespace numerical_methods;

// Create augmented matrix [A|b]
std::vector<std::vector<double>> augmented = {
    {2, 4, 2, 1, 10},
    {2, 2, 3, 3, 6},
    {4, 2, 2, 1, 6},
    {0, 2, 1, 1, 4}
};

// Create solver and solve
GaussianSolver solver(std::move(augmented));
auto solution = solver.solve();

// Verify solution
double residual = solver.verify_solution(coeff_matrix, rhs);
```

## Input Data Format

The program reads data from `data.txt` with the following format:

```
n
a₁₁ a₁₂ ... a₁ₙ b₁
a₂₁ a₂₂ ... a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```

### Example `data.txt`:
```
4
2 4 2 1 10
2 2 3 3 6
4 2 2 1 6
0 2 1 1 4
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
Gaussian Elimination Solver
============================================================
Successfully loaded 4x4 system from data.txt

Original System of Equations:
------------------------------------------------------------
     2.000*x1 +   4.000*x2 +   2.000*x3 +   1.000*x4 =     10.000
     2.000*x1 +   2.000*x2 +   3.000*x3 +   3.000*x4 =      6.000
     4.000*x1 +   2.000*x2 +   2.000*x3 +   1.000*x4 =      6.000
     0.000*x1 +   2.000*x2 +   1.000*x3 +   1.000*x4 =      4.000
------------------------------------------------------------

============================================================
GAUSSIAN ELIMINATION RESULTS
============================================================

System size: 4 equations, 4 unknowns

Solution:
------------------------------
x 1 =     -1.00000000
x 2 =      1.00000000
x 3 =      6.00000000
x 4 =     -4.00000000

Solution Verification:
------------------------------
Maximum residual: 1.77635684e-15
✓ Solution is ACCURATE

============================================================
```

## Algorithm Complexity

- **Time Complexity**: O(n³) for elimination + O(n²) for back substitution = O(n³)
- **Space Complexity**: O(n²) for the augmented matrix
- **Numerical Stability**: Enhanced through partial pivoting

## Error Handling

The implementation provides comprehensive error handling:

- **Matrix Validation**: Size consistency, empty matrix detection
- **Singular Matrix Detection**: Prevents division by zero
- **Input Validation**: File format checking, non-finite value detection
- **Numerical Issues**: Condition number monitoring

## Performance Optimizations

1. **Partial Pivoting**: Improves numerical stability
2. **Early Zero Detection**: Skips unnecessary operations when factors are effectively zero
3. **Memory Pre-allocation**: Reduces memory allocations with `reserve()`
4. **Vectorized Operations**: Uses `std::inner_product` for residual computation
5. **Const Correctness**: Enables compiler optimizations
6. **Cached Values**: Stores pivot values to avoid repeated calculations

## Mathematical Properties

- **Uniqueness**: For non-singular matrices, solution is unique
- **Stability**: Partial pivoting ensures numerical stability
- **Condition Number**: Measures sensitivity to input perturbations
- **Determinant**: Computed as product of diagonal elements after elimination

## Solution Verification

The implementation includes automatic solution verification:
- **Residual Computation**: Computes ||Ax - b||∞ for verification
- **Accuracy Assessment**: Categorizes solution accuracy based on residual magnitude
- **Visual Indicators**: Uses symbols (✓, ⚠, ✗) to indicate solution quality

## Condition Number Analysis

The solver provides condition number estimation:
- **Well-conditioned**: κ < 10¹² (reliable solution)
- **Ill-conditioned**: κ ≥ 10¹² (potentially unreliable)
- **Singular**: κ = ∞ (no unique solution)

## Applications

- **Engineering**: Structural analysis, circuit analysis, fluid dynamics
- **Scientific Computing**: Numerical simulations, data fitting
- **Computer Graphics**: Transformation matrices, lighting calculations
- **Economics**: Input-output models, linear programming
- **Machine Learning**: Linear regression, neural network training
- **Statistics**: Least squares fitting, regression analysis

## Limitations

- **Singular Matrices**: Cannot solve systems with zero determinant
- **Ill-conditioned Systems**: May produce inaccurate results for poorly conditioned matrices
- **Memory Usage**: O(n²) space requirement becomes significant for large systems
- **Computational Cost**: O(n³) time complexity limits scalability for very large systems

## Numerical Stability Considerations

1. **Partial Pivoting**: Reduces accumulation of round-off errors
2. **Epsilon Tolerance**: Uses machine epsilon scaled by safety factor
3. **Condition Monitoring**: Warns about potentially unstable computations
4. **Residual Checking**: Validates solution accuracy

## Alternative Methods

For specific use cases, consider:
- **LU Decomposition**: More efficient for multiple right-hand sides
- **QR Decomposition**: Better numerical stability for ill-conditioned systems
- **Iterative Methods**: More memory-efficient for large sparse systems
- **Specialized Solvers**: For matrices with special structure (tridiagonal, symmetric, etc.)

## Testing and Validation

The implementation includes comprehensive testing:
- **Known Solutions**: Verification against analytically solvable systems
- **Edge Cases**: Handling of singular and near-singular matrices
- **Numerical Precision**: Residual analysis for accuracy assessment
- **Performance**: Timing analysis for complexity verification

## Future Enhancements

- **LU Decomposition**: Implement LU factorization for improved efficiency
- **Iterative Refinement**: Add iterative improvement for better accuracy
- **Parallel Processing**: Multi-threaded elimination for large matrices
- **Sparse Matrix Support**: Memory-efficient handling of sparse systems
- **Advanced Pivoting**: Implement complete pivoting for maximum stability

