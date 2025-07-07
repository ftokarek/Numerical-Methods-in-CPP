# Jacobi Iterative Method Implementation

## Overview

This project implements the **Jacobi iterative method** in C++ using modern, professional techniques. 
The Jacobi method is an iterative algorithm for solving systems of linear equations Ax = b, particularly effective for large sparse systems where direct methods become computationally expensive.

## Mathematical Background

The Jacobi method solves the system Ax = b by decomposing the matrix A into:
```
A = D + L + U
```

Where:
- **D** is the diagonal matrix
- **L** is the strictly lower triangular matrix  
- **U** is the strictly upper triangular matrix

### Jacobi Iteration Formula

The iterative formula is:
```
x_i^(k+1) = (b_i - Σ(j≠i) a_ij * x_j^(k)) / a_ii
```

Or in matrix form:
```
x^(k+1) = D^(-1)(b - (L + U)x^(k))
```

### Convergence Conditions

The Jacobi method is **guaranteed to converge** if the matrix A is:
1. **Strictly diagonally dominant**: |a_ii| > Σ(j≠i)|a_ij| for all i
2. **Symmetric positive definite**
3. **Irreducibly diagonally dominant**

## Features

- **Convergence Analysis**: Automatic diagonal dominance checking
- **Adaptive Convergence**: Multiple convergence criteria (solution difference and residual)
- **Robust Error Handling**: Input validation and numerical stability checks
- **Performance Monitoring**: Iteration counting and residual tracking
- **Flexible Parameters**: Customizable tolerance and maximum iterations
- **Solution Verification**: Residual computation for accuracy assessment

## Files Structure

```
05_jacobi_method/
├── README.md                    # This file
├── main.cpp                     # Main program with file input
├── jacobi_solver.hpp            # Header with JacobiSolver class
├── data.txt                     # Sample input data file
└── main                         # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::JacobiSolver`
```cpp
class JacobiSolver 
{
    public:
        JacobiSolver(std::vector<std::vector<double>> matrix, std::vector<double> rhs);
        [[nodiscard]] std::vector<double> solve(
            const std::vector<double>& initial_guess = {},
            double tolerance = 1e-10,
            std::size_t max_iterations = 1000
        );
        [[nodiscard]] std::pair<bool, bool> analyze_convergence() const;
        [[nodiscard]] double verify_solution() const;
        [[nodiscard]] std::size_t get_iterations_performed() const noexcept;
        [[nodiscard]] double get_final_residual() const noexcept;
};
```

## Usage Example

```cpp
#include "jacobi_solver.hpp"
using namespace numerical_methods;

// Create system matrix and RHS
std::vector<std::vector<double>> matrix = {
    {10, -1,  2},
    {-1, 11, -1},
    { 2, -1, 10}
};
std::vector<double> rhs = {6, 25, -11};

// Create solver and solve
JacobiSolver solver(std::move(matrix), std::move(rhs));
auto solution = solver.solve();

// Check convergence and verify
auto [is_dominant, guaranteed] = solver.analyze_convergence();
double residual = solver.verify_solution();
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
10 -1 2
-1 11 -1
2 -1 10
6 25 -11
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
Jacobi Iterative Method Solver
============================================================
Successfully loaded 3x3 system from data.txt

System of equations (augmented matrix):
------------------------------------------------------------
   10.00    -1.00     2.00 |     6.00
   -1.00    11.00    -1.00 |    25.00
    2.00    -1.00    10.00 |   -11.00
------------------------------------------------------------

Convergence Analysis:
------------------------------
✓ Matrix is diagonally dominant
✓ Convergence is GUARANTEED

Use default parameters (tolerance=1e-10, max_iter=1000)? (y/n): y

Solving system with Jacobi method...

============================================================
JACOBI METHOD RESULTS
============================================================

System size: 3x3
Iterations performed: 19

Solution:
------------------------------
x 1 =      1.04326923
x 2 =      2.26923077
x 3 =     -1.08173077

Solution Verification:
------------------------------
Final residual: 1.18854260e-10
⚠ Solution has acceptable accuracy

============================================================
```

## Algorithm Complexity

- **Time Complexity per Iteration**: O(n²)
- **Space Complexity**: O(n²) for matrix storage
- **Convergence Rate**: Linear (geometric progression)
- **Memory Access**: Cache-friendly for sparse matrices

## Convergence Analysis

### Diagonal Dominance Types
1. **Strictly Diagonally Dominant**: Convergence guaranteed
2. **Weakly Diagonally Dominant**: Convergence possible
3. **Not Diagonally Dominant**: Convergence unlikely

### Convergence Criteria
The method stops when either:
1. **Solution difference**: ||x^(k+1) - x^(k)||∞ < tolerance
2. **Maximum iterations**: Reached iteration limit

## Advantages and Disadvantages

### Advantages
- **Parallelizable**: Each component can be computed independently
- **Memory Efficient**: Only requires matrix storage, no additional factorization
- **Numerically Stable**: No pivoting or matrix modifications needed
- **Simple Implementation**: Straightforward iterative formula

### Disadvantages
- **Slow Convergence**: Linear convergence rate
- **Convergence Requirements**: Needs diagonally dominant matrices
- **Not Universal**: Cannot solve all linear systems
- **Iteration Dependent**: Number of iterations varies with condition

## Applications

### Scientific Computing
- **Finite Difference Methods**: Discretized PDEs
- **Image Processing**: Iterative filtering and restoration
- **Computational Fluid Dynamics**: Pressure correction algorithms

### Engineering
- **Structural Analysis**: Large sparse systems
- **Electrical Networks**: Circuit analysis
- **Heat Transfer**: Steady-state temperature distribution

### Specialized Domains
- **Graph Theory**: Network flow problems
- **Machine Learning**: Distributed optimization
- **Numerical PDEs**: Relaxation methods

## Error Handling

The implementation provides comprehensive error handling:

- **Matrix Validation**: Square matrix, diagonal non-zero elements
- **Convergence Analysis**: Diagonal dominance checking
- **Input Validation**: File format and numerical validity
- **Iteration Monitoring**: Maximum iteration limits and divergence detection

## Performance Optimizations

1. **Early Termination**: Exact node matching for interpolation points
2. **Vectorized Operations**: Use of `std::inner_product` for residual computation
3. **Memory Pre-allocation**: Reserved vectors to avoid reallocation
4. **Const Correctness**: Enables compiler optimizations
5. **Move Semantics**: Efficient resource management

## Mathematical Properties

- **Convergence Rate**: O(ρ^k) where ρ is spectral radius
- **Stability**: Unconditionally stable for diagonally dominant matrices
- **Accuracy**: Limited by tolerance and matrix condition number
- **Robustness**: Handles ill-conditioned systems better than direct methods
