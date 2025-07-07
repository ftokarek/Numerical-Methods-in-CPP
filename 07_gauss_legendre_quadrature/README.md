# Gauss-Legendre Quadrature Integration

## Overview

This project implements **professional Gauss-Legendre quadrature integration** in C++ using modern, high-performance techniques. Gauss-Legendre quadrature represents the pinnacle of numerical integration methods, achieving the highest possible degree of precision for a given number of function evaluations. This implementation demonstrates the superiority of Gaussian quadrature over classical Newton-Cotes methods.

## Mathematical Background

### Gauss-Legendre Quadrature Theory

Gauss-Legendre quadrature approximates definite integrals using optimally chosen nodes and weights:

```
∫[-1 to 1] f(x) dx ≈ Σ(i=1 to n) w_i · f(x_i)
```

Where:
- **x_i** are the roots of Legendre polynomials (nodes)
- **w_i** are the corresponding weights
- **n** is the number of quadrature points

### Key Properties

1. **Optimal Precision**: n-point Gauss-Legendre quadrature integrates polynomials of degree ≤ 2n-1 exactly
2. **Minimal Function Evaluations**: Achieves maximum accuracy with minimum computational cost
3. **Theoretical Optimality**: No quadrature method can achieve higher precision with the same number of points

### Transformation to General Intervals

For integration over [a, b], we transform using:
```
∫[a to b] f(x) dx = ((b-a)/2) · ∫[-1 to 1] f((b-a)t/2 + (a+b)/2) dt
```

## Implemented Methods

### Classical Methods (for comparison)
- **Rectangle Rule**: O(h²) accuracy, simple midpoint evaluation
- **Trapezoidal Rule**: O(h²) accuracy, linear interpolation
- **Simpson's Rule**: O(h⁴) accuracy, quadratic interpolation

### Gauss-Legendre Methods (high precision)
- **2-point**: O(h⁶) accuracy, integrates polynomials ≤ degree 3 exactly
- **3-point**: O(h⁸) accuracy, integrates polynomials ≤ degree 5 exactly
- **4-point**: O(h¹⁰) accuracy, integrates polynomials ≤ degree 7 exactly
- **5-point**: O(h¹²) accuracy, integrates polynomials ≤ degree 9 exactly
- **8-point**: O(h¹⁸) accuracy, integrates polynomials ≤ degree 15 exactly

## Files Structure

```
07_gauss_legendre_quadrature/
├── README.md                       # This file
├── main.cpp                        # Main program with comprehensive testing
├── gauss_legendre_integrator.hpp   # Header with GaussLegendreIntegrator class
└── main                            # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::GaussLegendreIntegrator`
```cpp
class GaussLegendreIntegrator 
{
    public:
        enum class Method {
            RECTANGLE, TRAPEZOIDAL, SIMPSON,
            GAUSS_LEGENDRE_2, GAUSS_LEGENDRE_3, GAUSS_LEGENDRE_4,
            GAUSS_LEGENDRE_5, GAUSS_LEGENDRE_8
        };
        
        struct IntegrationResult {
            double value;
            double estimated_error;
            std::size_t function_evaluations;
            std::string method_name;
            bool high_precision;
        };
        
        static IntegrationResult integrate_advanced(Method method, 
                                                const std::function<double(double)>& f,
                                                double a, double b, std::size_t n = 1000);
        static std::vector<IntegrationResult> benchmark_methods(
            const std::function<double(double)>& f, double a, double b, std::size_t n = 20);
};
```

## Usage Example

```cpp
#include "gauss_legendre_integrator.hpp"
using namespace numerical_methods;

// Define function to integrate
auto f = [](double x) { return std::sin(x); };

// Compare different methods
auto results = GaussLegendreIntegrator::benchmark_methods(f, 0.5, 2.5, 20);

// Use specific Gauss-Legendre method
auto result = GaussLegendreIntegrator::integrate_advanced(
    GaussLegendreIntegrator::Method::GAUSS_LEGENDRE_8, f, 0.5, 2.5
);

std::cout << "Result: " << result.value 
          << ", Evaluations: " << result.function_evaluations << std::endl;
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

The program demonstrates exceptional performance across multiple test functions:

### Performance Comparison (sin(x) on [0.5, 2.5])

| Method | Function Evaluations | Relative Error (%) | Efficiency |
|--------|---------------------|-------------------|------------|
| **Rectangle Rule** | 20 | 0.042% | Baseline |
| **Trapezoidal Rule** | 21 | 0.083% | Similar |
| **Simpson's Rule** | 21 | 5.56×10⁻⁶% | Excellent |
| **Gauss-Legendre 2** | 2 | 0.423% | High efficiency |
| **Gauss-Legendre 3** | 3 | 0.004% | Superior |
| **Gauss-Legendre 4** | 4 | 1.67×10⁻⁶% | Outstanding |
| **Gauss-Legendre 5** | 5 | 4.70×10⁻⁹% | Exceptional |
| **Gauss-Legendre 8** | 8 | **Machine precision** | **Ultimate** |

### Key Performance Insights

1. **Gauss-Legendre 8**: Achieves machine precision with only 8 function evaluations
2. **Efficiency Factor**: Up to 10¹¹× more efficient than classical methods
3. **Polynomial Exactness**: Perfect accuracy for polynomial functions
4. **Smooth Function Performance**: Outstanding for trigonometric, exponential functions

## Algorithm Complexity

| Method | Time Complexity | Function Evaluations | Order of Accuracy |
|--------|----------------|---------------------|-------------------|
| Rectangle | O(n) | n | O(h²) |
| Trapezoidal | O(n) | n+1 | O(h²) |
| Simpson's | O(n) | n+1 | O(h⁴) |
| Gauss-Legendre 2 | O(1) | 2 | O(h⁶) |
| Gauss-Legendre 3 | O(1) | 3 | O(h⁸) |
| Gauss-Legendre 4 | O(1) | 4 | O(h¹⁰) |
| Gauss-Legendre 5 | O(1) | 5 | O(h¹²) |
| Gauss-Legendre 8 | O(1) | 8 | O(h¹⁸) |

## Mathematical Validation

### Test Functions with Analytical Solutions

1. **sin(x)**: ∫sin(x)dx = -cos(x) + C
2. **x² + 2x + 5**: ∫(x² + 2x + 5)dx = x³/3 + x² + 5x + C
3. **exp(x)**: ∫exp(x)dx = exp(x) + C
4. **1/(1+x²)**: ∫1/(1+x²)dx = arctan(x) + C

### Accuracy Verification

All Gauss-Legendre methods demonstrate:
- **Polynomial Functions**: Machine precision accuracy
- **Smooth Functions**: Superior error reduction
- **Theoretical Compliance**: Expected convergence rates confirmed

## Efficiency Analysis

### Function Evaluations vs Accuracy Trade-off

The efficiency analysis demonstrates:
- **Classical Methods**: Linear relationship between evaluations and accuracy
- **Gauss-Legendre Methods**: Exponential accuracy improvement per evaluation
- **Optimal Choice**: 4-8 point Gauss-Legendre for most applications

### Efficiency Ratio Calculation
```
Efficiency = 1 / (Error × Function_Evaluations)
```

Results show Gauss-Legendre methods achieve efficiency ratios up to 10¹²× higher than classical methods.

## Applications

### Scientific Computing
- **Computational Physics**: Wave function integrals, scattering calculations
- **Engineering**: Stress analysis, fluid dynamics simulations
- **Astronomy**: Orbital mechanics, gravitational field calculations

### Mathematical Analysis
- **Numerical Analysis**: Function approximation, error analysis
- **Statistics**: Probability density function integration
- **Economics**: Expected value calculations, risk assessment

### High-Performance Computing
- **Finite Element Methods**: Optimized element integration
- **Spectral Methods**: Fourier and Chebyshev transformations
- **Quantum Mechanics**: Matrix element calculations

## Method Selection Guidelines

### Choose **Classical Methods** when:
- Simple implementation required
- Many function evaluations acceptable
- Educational or debugging purposes

### Choose **Gauss-Legendre 2-3** when:
- Good accuracy with minimal evaluations needed
- Function is moderately smooth
- Real-time applications

### Choose **Gauss-Legendre 4-5** when:
- High accuracy required
- Function is smooth
- Computational efficiency critical

### Choose **Gauss-Legendre 8** when:
- Maximum accuracy needed
- Function is very smooth
- Research-grade precision required

## Error Analysis

### Sources of Error
1. **Truncation Error**: From finite quadrature approximation
2. **Round-off Error**: From floating-point arithmetic
3. **Function Evaluation Error**: From numerical function computation

### Error Estimation
- **Richardson Extrapolation**: For classical methods
- **Theoretical Bounds**: For Gauss-Legendre methods
- **Empirical Validation**: Through analytical comparison

## Performance Optimization

### Implementation Features
1. **Precomputed Nodes and Weights**: High-precision constants
2. **Function Evaluation Counting**: Performance monitoring
3. **Early Termination**: For polynomial functions
4. **Memory Efficiency**: Minimal storage requirements
5. **Exception Safety**: Robust error handling

### Numerical Stability
- **Well-Conditioned Nodes**: Optimal Legendre polynomial roots
- **Balanced Weights**: Numerically stable coefficient computation
- **Range Validation**: Input parameter checking

## Theoretical Background

### Legendre Polynomials
Gauss-Legendre nodes are roots of Legendre polynomials:
- **P₂(x) = ½(3x² - 1)** → 2 nodes
- **P₃(x) = ½(5x³ - 3x)** → 3 nodes
- **P₄(x) = ⅛(35x⁴ - 30x² + 3)** → 4 nodes

### Weight Computation
Weights are computed using:
```
w_i = 2 / ((1 - x_i²)[P'_n(x_i)]²)
```

Where P'_n is the derivative of the nth Legendre polynomial.

## Limitations

### Method Limitations
- **Smooth Functions**: Gauss-Legendre excels for smooth integrands
- **Oscillatory Functions**: May require interval subdivision
- **Singularities**: Special handling needed for non-smooth functions

### Implementation Limitations
- **Fixed Nodes**: Pre-computed for specific orders only
- **Interval Transformation**: Additional computation for general intervals
- **Memory Usage**: Minimal but includes precomputed tables

