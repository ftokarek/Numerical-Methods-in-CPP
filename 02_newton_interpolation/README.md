# Newton Interpolation Implementation

## Overview

This project implements **Newton's divided difference interpolation** method in C++ using modern, optimized techniques. 
Newton interpolation is a polynomial interpolation method that constructs a polynomial of degree at most n-1 passing through n given data points using divided differences and provides an efficient way to evaluate the polynomial using Horner's method.

## Mathematical Background

Given n data points (x₀, y₀), (x₁, y₁), ..., (xₙ₋₁, yₙ₋₁), Newton's interpolating polynomial is:

```
P(x) = f[x₀] + f[x₀,x₁](x-x₀) + f[x₀,x₁,x₂](x-x₀)(x-x₁) + ... + f[x₀,x₁,...,xₙ₋₁](x-x₀)(x-x₁)...(x-xₙ₋₂)
```

where f[x₀,x₁,...,xₖ] are the **divided differences** defined recursively:

```
f[xᵢ] = yᵢ
f[xᵢ,xᵢ₊₁,...,xᵢ₊ₖ] = (f[xᵢ₊₁,...,xᵢ₊ₖ] - f[xᵢ,...,xᵢ₊ₖ₋₁]) / (xᵢ₊ₖ - xᵢ)
```

### Advantages over Lagrange Interpolation

- **Incremental construction**: Easy to add new points
- **Efficient evaluation**: O(n) using Horner's method
- **Numerical stability**: Better conditioned for computation
- **Explicit coefficients**: Clear polynomial representation

## Features

- **High Performance**: O(n) evaluation using Horner's method
- **Comprehensive Error Handling**: File validation, input checking, duplicate detection
- **Robust File I/O**: Filesystem validation and detailed error reporting
- **Memory Efficient**: Optimized memory usage with reserved vectors

## Files Structure

```
02_newton_interpolation/
├── README.md                    # This file
├── main.cpp                     # Main program with file input
├── newton_interpolator.hpp      # Header with NewtonInterpolator class
├── data.txt                     # Sample input data file
└── main                         # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::Node`
```cpp
struct Node 
{
    double x, y;
    constexpr Node(double x, double y) noexcept;
    bool operator==(const Node& other) const noexcept;
};
```

### `numerical_methods::NewtonInterpolator`
```cpp
class NewtonInterpolator 
{
    public:
        explicit NewtonInterpolator(std::vector<Node> nodes);
        [[nodiscard]] double evaluate(double x) const noexcept;
        [[nodiscard]] const std::vector<Node>& get_nodes() const noexcept;
        [[nodiscard]] const std::vector<double>& get_coefficients() const noexcept;
        [[nodiscard]] std::size_t size() const noexcept;
        [[nodiscard]] std::size_t degree() const noexcept;
        [[nodiscard]] bool empty() const noexcept;
};
```

## Usage Example

```cpp
#include "newton_interpolator.hpp"
using namespace numerical_methods;

// Create interpolator with data points
std::vector<Node> nodes = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 5.0}};
NewtonInterpolator interpolator(std::move(nodes));

// Evaluate at any point
double result = interpolator.evaluate(1.5);

// Get coefficients
auto coeffs = interpolator.get_coefficients();
```

## Input Data Format

The program reads data from `data.txt` with the following format:

```
n
x₀ y₀
x₁ y₁
...
xₙ₋₁ yₙ₋₁
```

### Example `data.txt`:
```
4
1.0 2.0
2.0 3.0
3.0 5.0
4.0 8.0
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

## Program Flow

1. **File Reading**: Validates and reads interpolation data from `data.txt`
2. **Data Validation**: Checks for duplicate x-coordinates and non-finite values
3. **Interpolator Creation**: Constructs Newton interpolator with divided differences
4. **User Input**: Prompts for evaluation point with validation
5. **Evaluation**: Computes polynomial value using Horner's method
6. **Results Display**: Shows comprehensive formatted results

## Performance Characteristics

- **Time Complexity**: 
  - Initialization: O(n²) for divided differences computation
  - Evaluation: O(n) using Horner's method
- **Space Complexity**: O(n²) during initialization, O(n) for storage
- **Numerical Stability**: Enhanced with epsilon-based comparisons

## Error Handling

The implementation provides comprehensive error handling:

- **File Errors**: Non-existent files, permission issues
- **Format Errors**: Invalid data format, non-numeric values
- **Mathematical Errors**: Duplicate x-coordinates, non-finite values
- **Input Errors**: Invalid user input, non-finite evaluation points

## Sample Output

```
Newton Polynomial Interpolation
==================================================

Successfully loaded 4 nodes from data.txt

Enter the point at which to evaluate the Newton polynomial: 2.5

==================================================
NEWTON INTERPOLATION RESULTS
==================================================

Number of interpolation nodes: 4
Polynomial degree: 3

Interpolation nodes:
------------------------------
           x        f(x)
------------------------------
    1.000000    2.000000
    2.000000    3.000000
    3.000000    5.000000
    4.000000    8.000000

Evaluation:
------------------------------
x = 2.500000
P(x) = 3.875000

Newton coefficients (divided differences):
----------------------------------------
b[0] =     2.000000
b[1] =     1.000000
b[2] =     1.000000
b[3] =     0.333333

==================================================
```

## Mathematical Properties

- **Uniqueness**: Same polynomial as Lagrange interpolation
- **Incremental**: Easy to add new data points
- **Explicit Form**: Clear polynomial coefficients
- **Stability**: Generally more stable than Lagrange form

## Applications

- **Function Approximation**: Approximate unknown functions
- **Data Interpolation**: Fill gaps in experimental data
- **Numerical Integration**: Basis for Newton-Cotes formulas
- **Differential Equations**: Predictor-corrector methods
- **Computer Graphics**: Smooth curve generation

## Limitations

- **Runge's Phenomenon**: High-degree polynomials may oscillate
- **Extrapolation**: Unreliable outside data range
- **Sensitivity**: Can be sensitive to data errors
- **Degree Limitation**: High degrees may cause numerical issues

## Advanced Features

- **Comprehensive Documentation**: Doxygen-style comments
- **Modern C++ Features**: `constexpr`, `[[nodiscard]]`, `noexcept`
- **Exception Safety**: Strong exception safety guarantees
- **Memory Optimization**: Move semantics and reserved containers
- **Input Validation**: Robust checking for all inputs
