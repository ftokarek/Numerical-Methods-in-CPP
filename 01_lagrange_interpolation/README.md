# Lagrange Interpolation Implementation

## Overview

This code implements the **Lagrange interpolation** method in C++ using both traditional and optimized barycentric forms. 
Lagrange interpolation is a polynomial interpolation technique that finds a polynomial of degree at most n-1 passing through n given data points.

## Mathematical Background

Given n data points (x₀, y₀), (x₁, y₁), ..., (xₙ₋₁, yₙ₋₁), the Lagrange interpolating polynomial is:

```
P(x) = Σ(i=0 to n-1) yᵢ * Lᵢ(x)
```

where Lᵢ(x) is the i-th Lagrange basis polynomial:

```
Lᵢ(x) = Π(j=0 to n-1, j≠i) (x - xⱼ)/(xᵢ - xⱼ)
```

### Barycentric Form

This implementation uses the **barycentric form** for improved efficiency:

```
P(x) = Σ(i=0 to n-1) wᵢ/(x - xᵢ) * yᵢ / Σ(i=0 to n-1) wᵢ/(x - xᵢ)
```

where wᵢ are the barycentric weights computed once during initialization.

## Features

- **High Performance**: O(n) evaluation complexity using barycentric form
- **Numerical Stability**: Epsilon-based floating-point comparisons
- **Professional C++ Design**: 
  - Class-based architecture
  - Exception handling
  - Const correctness
  - Move semantics
  - Modern C++ attributes (`[[nodiscard]]`, `noexcept`)
- **Memory Efficient**: Pre-computed weights and reserved vectors
- **Early Termination**: Direct return for exact node matches

## Files Structure

```
01_lagrange_interpolation/
├── README.md                    # This file
├── main.cpp                     # Main program with two demonstration tasks
├── lagrange_interpolator.hpp    # Header file with LagrangeInterpolator class
└── main                         # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::Node`
```cpp
struct Node {
    double x, y;
    constexpr Node(double x, double y) noexcept;
};
```

### `numerical_methods::LagrangeInterpolator`
```cpp
class LagrangeInterpolator {
public:
    explicit LagrangeInterpolator(std::vector<Node> nodes);
    [[nodiscard]] double evaluate(double x) const noexcept;
    [[nodiscard]] const std::vector<Node>& get_nodes() const noexcept;
    [[nodiscard]] size_t size() const noexcept;
};
```

## Usage Example

```cpp
#include "lagrange_interpolator.hpp"
using namespace numerical_methods;

// Create interpolator with data points
LagrangeInterpolator interpolator({{-4, 5}, {-3, 2}, {1, 5}, {2, 2}});

// Evaluate at any point
double result = interpolator.evaluate(0.5);
```

## Compilation

```bash
g++ -std=c++17 -O3 -Wall -Wextra main.cpp -o main
```
# or

```bash
clang++ -std=c++17 main.cpp -o  main
```

```bash
# Run the program
./main
```

## Program Tasks

The main program demonstrates two interpolation tasks:

### Task 1: Interactive Interpolation
- Uses nodes: (-4, 5), (-3, 2), (1, 5), (2, 2)
- Prompts user for evaluation point
- Displays interpolation result

### Task 2: Cube Root Approximation
- Uses nodes: (27, 3), (64, 4), (125, 5), (216, 6)
- Approximates cube root of 50
- Expected result: ≈ 3.684

## Performance Characteristics

- **Time Complexity**: 
  - Initialization: O(n²)
  - Evaluation: O(n)
- **Space Complexity**: O(n)
- **Numerical Stability**: Enhanced with epsilon-based comparisons

## Error Handling

The implementation includes robust error handling:
- Throws `std::invalid_argument` for empty node sets
- Throws `std::invalid_argument` for duplicate x-coordinates
- Uses epsilon tolerance for floating-point comparisons

## Mathematical Properties

- **Exactness**: The polynomial passes exactly through all given points
- **Uniqueness**: There exists exactly one polynomial of degree ≤ n-1 through n points
- **Oscillation**: May exhibit Runge's phenomenon for high-degree polynomials

## Limitations

- **Runge's Phenomenon**: High-degree polynomials may oscillate wildly
- **Extrapolation**: Results outside the data range may be unreliable
- **Condition Number**: Poorly conditioned for closely spaced points

## Applications

- Function approximation
- Data smoothing
- Numerical integration (quadrature rules)
- Computer graphics (curve fitting)
- Signal processing

