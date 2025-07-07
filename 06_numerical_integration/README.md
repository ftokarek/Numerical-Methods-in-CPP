# Numerical Integration Implementation

## Overview

This project implements **professional numerical integration methods** in C++ using modern, optimized techniques. The implementation provides multiple quadrature methods including Rectangle (Midpoint), Trapezoidal, Simpson's 1/3 rule, and Adaptive Simpson integration with comprehensive error analysis and performance monitoring.

## Mathematical Background

Numerical integration approximates the definite integral:

```
∫[a to b] f(x) dx ≈ Σ w_i · f(x_i)
```

Where w_i are weights and x_i are evaluation points determined by the chosen method.

### Implemented Methods

#### 1. Rectangle (Midpoint) Rule - O(h²)
```
∫[a to b] f(x) dx ≈ h · Σ f(x_i + h/2)
```
- Uses function values at interval midpoints
- Second-order accuracy
- Excellent for smooth functions

#### 2. Trapezoidal Rule - O(h²)
```
∫[a to b] f(x) dx ≈ h/2 · [f(a) + 2·Σf(x_i) + f(b)]
```
- Linear interpolation between points
- Second-order accuracy
- Good general-purpose method

#### 3. Simpson's 1/3 Rule - O(h⁴)
```
∫[a to b] f(x) dx ≈ h/3 · [f(a) + 4·Σf(x_odd) + 2·Σf(x_even) + f(b)]
```
- Quadratic interpolation
- Fourth-order accuracy
- Superior for polynomial functions

#### 4. Adaptive Simpson Rule
- Automatic subdivision based on error estimates
- Richardson extrapolation for error control
- Machine precision accuracy for smooth functions

## Features

- **Multiple Integration Methods**: Rectangle, Trapezoidal, Simpson's, Adaptive Simpson
- **Comprehensive Error Analysis**: Richardson extrapolation and convergence monitoring
- **Performance Optimization**: Efficient algorithms with minimal function evaluations
- **Validation Suite**: Test functions with known analytical solutions
- **Convergence Analysis**: Demonstrates theoretical convergence rates
- **Interactive Testing**: User-customizable parameters

## Files Structure

```
06_numerical_integration/
├── README.md                    # This file
├── main.cpp                     # Main program with comprehensive testing
├── numerical_integrator.hpp     # Header with NumericalIntegrator class
└── main                         # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::NumericalIntegrator`
```cpp
class NumericalIntegrator 
{
    public:
        enum class Method 
        {
            RECTANGLE, TRAPEZOIDAL, SIMPSON, ADAPTIVE_SIMPSON
        };
        
        static double integrate(Method method, const std::function<double(double)>& f,
                            double a, double b, std::size_t n = 1000,
                            double tolerance = 1e-8);
        static double estimate_error(Method method, const std::function<double(double)>& f,
                                    double a, double b, std::size_t n);
        static std::size_t get_order_of_accuracy(Method method);
        static std::string get_method_name(Method method);
};
```

## Usage Example

```cpp
#include "numerical_integrator.hpp"
using namespace numerical_methods;

// Define function to integrate
auto f = [](double x) { return std::sin(x); };

// Integrate using different methods
double result_simpson = NumericalIntegrator::integrate(
    NumericalIntegrator::Method::SIMPSON, f, 0.5, 2.5, 100
);

double result_adaptive = NumericalIntegrator::integrate(
    NumericalIntegrator::Method::ADAPTIVE_SIMPSON, f, 0.5, 2.5, 0, 1e-10
);

// Estimate integration error
double error = NumericalIntegrator::estimate_error(
    NumericalIntegrator::Method::SIMPSON, f, 0.5, 2.5, 100
);
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

## Sample Output Analysis

The program demonstrates excellent accuracy across multiple test functions:

### Test Function Results (100 intervals, [0.5, 2.5])

| Function | Method | Relative Error (%) | Performance |
|----------|--------|--------------------|-------------|
| **sin(x)** | Rectangle | 0.0017% | Good |
| | Trapezoidal | 0.0033% | Good |
| | Simpson's | 8.89e-6% | Excellent |
| | Adaptive | Machine precision | Perfect |
| **x² + 2x + 5** | Rectangle | 0.00031% | Excellent |
| | Trapezoidal | 0.00063% | Excellent |
| | Simpson's | Machine precision | Perfect |
| | Adaptive | Machine precision | Perfect |
| **exp(x)** | Rectangle | 0.0017% | Good |
| | Trapezoidal | 0.0033% | Good |
| | Simpson's | 8.89e-6% | Excellent |
| | Adaptive | Machine precision | Perfect |
| **1/(1+x²)** | Rectangle | 0.00125% | Excellent |
| | Trapezoidal | 0.0025% | Good |
| | Simpson's | 4.65e-5% | Excellent |
| | Adaptive | Machine precision | Perfect |

### Convergence Analysis Verification

The convergence analysis confirms theoretical expectations:

- **Rectangle & Trapezoidal**: O(h²) - errors decrease by ~4× when intervals double
- **Simpson's Rule**: O(h⁴) - errors decrease by ~16× when intervals double

## Algorithm Complexity

| Method | Time Complexity | Space Complexity | Order of Accuracy |
|--------|----------------|------------------|-------------------|
| Rectangle | O(n) | O(1) | O(h²) |
| Trapezoidal | O(n) | O(1) | O(h²) |
| Simpson's | O(n) | O(1) | O(h⁴) |
| Adaptive | O(n log n) | O(log n) | O(h⁴) adaptive |

## Error Analysis

### Richardson Extrapolation
For methods with known convergence order p:
```
Error ≈ |I_{2n} - I_n| / (2^p - 1)
```

### Adaptive Error Control
The adaptive Simpson method uses recursive subdivision with error estimates:
```
Error ≈ |S_split - S_whole| / 15
```

## Test Functions with Analytical Solutions

1. **sin(x)**: ∫sin(x)dx = -cos(x) + C
2. **x² + 2x + 5**: ∫(x² + 2x + 5)dx = x³/3 + x² + 5x + C
3. **exp(x)**: ∫exp(x)dx = exp(x) + C
4. **1/(1+x²)**: ∫1/(1+x²)dx = arctan(x) + C

## Performance Characteristics

### Accuracy Hierarchy (for smooth functions):
1. **Adaptive Simpson**: Machine precision
2. **Simpson's Rule**: Excellent (10⁻⁶ to 10⁻¹² relative error)
3. **Rectangle Rule**: Good (10⁻³ to 10⁻⁴ relative error)
4. **Trapezoidal Rule**: Good (10⁻³ to 10⁻⁴ relative error)

### Computational Efficiency:
- **Function Evaluations**: Rectangle = Trapezoidal < Simpson's < Adaptive
- **Setup Cost**: All methods have O(1) setup
- **Memory Usage**: All methods use O(1) additional memory

## Applications

### Scientific Computing
- **Physics Simulations**: Energy calculations, flux integrals
- **Engineering**: Area under curves, work calculations
- **Statistics**: Probability density function integration

### Numerical Analysis
- **Differential Equations**: Numerical solution methods
- **Approximation Theory**: Function approximation
- **Computational Mathematics**: Algorithm validation

### Real-World Applications
- **Financial Mathematics**: Option pricing, risk assessment
- **Computer Graphics**: Lighting calculations, ray tracing
- **Signal Processing**: Fourier transforms, filtering

## Method Selection Guidelines

### Choose **Rectangle Rule** when:
- Simple implementation needed
- Function is smooth
- Moderate accuracy acceptable

### Choose **Trapezoidal Rule** when:
- Linear approximation is appropriate
- Endpoints are critical
- Legacy compatibility needed

### Choose **Simpson's Rule** when:
- High accuracy required
- Function is polynomial-like
- Fixed number of intervals

### Choose **Adaptive Simpson** when:
- Maximum accuracy needed
- Function smoothness varies
- Computational cost is secondary

## Error Handling

The implementation provides comprehensive error handling:

- **Parameter Validation**: Finite bounds, positive intervals
- **Function Evaluation**: Handles mathematical exceptions
- **Convergence Monitoring**: Adaptive method convergence checking
- **Numerical Stability**: Epsilon-based comparisons

## Limitations

### General Limitations:
- **Singular Functions**: Methods may fail near singularities
- **Highly Oscillatory Functions**: May require many intervals
- **Discontinuous Functions**: Reduced accuracy at discontinuities

### Method-Specific Limitations:
- **Rectangle/Trapezoidal**: Slow convergence for rough functions
- **Simpson's**: Requires even number of intervals
- **Adaptive**: May be computationally expensive

## Validation and Testing

The implementation includes comprehensive validation:
- **Analytical Verification**: Known exact solutions
- **Convergence Testing**: Theoretical rate confirmation
- **Error Estimation**: Richardson extrapolation validation
- **Edge Case Handling**: Boundary condition testing
