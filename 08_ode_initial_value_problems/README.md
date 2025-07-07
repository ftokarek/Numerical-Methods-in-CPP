# ODE Initial Value Problems Solver

## Overview

This project implements **professional ODE solvers** for initial value problems (Cauchy problems) in C++ using modern, high-performance techniques. The implementation provides multiple numerical methods for solving first-order ordinary differential equations of the form dy/dx = f(x, y) with initial condition y(x₀) = y₀.

## Mathematical Background

### Initial Value Problems (Cauchy Problems)

An initial value problem consists of:
```
dy/dx = f(x, y)
y(x₀) = y₀
```

Where:
- **f(x, y)** is a given function
- **x₀, y₀** are initial conditions
- **y(x)** is the unknown solution function

### Numerical Methods Implemented

#### 1. Euler's Method - O(h¹)
```
y_{n+1} = y_n + h·f(x_n, y_n)
```
- Simple first-order method
- Linear convergence rate
- Good for educational purposes

#### 2. Modified Euler (Heun's Method) - O(h²)
```
k₁ = f(x_n, y_n)
k₂ = f(x_n + h, y_n + h·k₁)
y_{n+1} = y_n + h·(k₁ + k₂)/2
```
- Second-order accuracy
- Predictor-corrector approach
- Better stability than Euler

#### 3. 2nd Order Runge-Kutta - O(h²)
```
k₁ = f(x_n, y_n)
k₂ = f(x_n + h/2, y_n + h·k₁/2)
y_{n+1} = y_n + h·k₂
```
- Midpoint method
- Second-order accuracy
- Good balance of accuracy and efficiency

#### 4. 4th Order Runge-Kutta - O(h⁴)
```
k₁ = f(x_n, y_n)
k₂ = f(x_n + h/2, y_n + h·k₁/2)
k₃ = f(x_n + h/2, y_n + h·k₂/2)
k₄ = f(x_n + h, y_n + h·k₃)
y_{n+1} = y_n + h·(k₁ + 2k₂ + 2k₃ + k₄)/6
```
- Fourth-order accuracy
- Industry standard for ODE solving
- Excellent accuracy-to-cost ratio

## Files Structure

```
08_ode_initial_value_problems/
├── README.md                    # This file
├── main.cpp                     # Main program with comprehensive testing
├── ode_solver.hpp               # Header with ODESolver class
└── ode                          # Compiled binary (ignored by git)
```

## Class Interface

### `numerical_methods::ODESolver`
```cpp
class ODESolver {
public:
    enum class Method {
        EULER, MODIFIED_EULER, RUNGE_KUTTA_2, RUNGE_KUTTA_4, ADAPTIVE_RK4
    };
    
    struct SolutionResult {
        std::vector<SolutionPoint> points;
        double final_value;
        std::size_t steps_taken;
        double estimated_error;
        std::string method_name;
        bool converged;
    };
    
    static double solve(Method method, const std::function<double(double, double)>& f,
                       double x0, double y0, double x_end, double h);
    static SolutionResult solve_complete(Method method, 
                                        const std::function<double(double, double)>& f,
                                        double x0, double y0, double x_end, double h);
};
```

## Usage Example

```cpp
#include "ode_solver.hpp"
using namespace numerical_methods;

// Define ODE: y' = 2y(x+1)
auto f = [](double x, double y) { return 2.0 * y * (x + 1.0); };

// Solve using RK4
double result = ODESolver::solve(ODESolver::Method::RUNGE_KUTTA_4, 
                                f, 0.0, 1.0, 1.0, 0.1);

// Get complete solution trajectory
auto solution = ODESolver::solve_complete(ODESolver::Method::RUNGE_KUTTA_4,
                                         f, 0.0, 1.0, 1.0, 0.1);
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

The program demonstrates exceptional performance across multiple test problems:

### Test Problem Results (h = 0.1)

#### Problem 1: y' = 2y(x+1), y(0) = 1

| Method | Result | Relative Error (%) | Performance |
|--------|--------|--------------------|-------------|
| **Euler** | 12.635 | 37.09% | Poor accuracy |
| **Modified Euler** | 19.306 | 3.88% | Good improvement |
| **RK2** | 19.101 | 4.90% | Similar to Modified Euler |
| **RK4** | 20.081 | **0.021%** | **Excellent accuracy** |

**Analytical**: y(1) = e³ ≈ 20.0855

#### Problem 2: y' = x + y, y(0) = 0.1

| Method | Result | Relative Error (%) | Performance |
|--------|--------|--------------------|-------------|
| **Euler** | 0.8531 | 13.84% | Moderate accuracy |
| **Modified Euler** | 0.9855 | 0.467% | Good accuracy |
| **RK2** | 0.9855 | 0.467% | Same as Modified Euler |
| **RK4** | 0.9901 | **0.0002%** | **Outstanding accuracy** |

**Analytical**: y(1) = -1 - 1 + 1.1e ≈ 0.9901

#### Problem 3: y' = -y, y(0) = 1

| Method | Result | Relative Error (%) | Performance |
|--------|--------|--------------------|-------------|
| **Euler** | 0.1216 | 10.17% | Acceptable |
| **Modified Euler** | 0.1358 | 0.36% | Good accuracy |
| **RK2** | 0.1358 | 0.36% | Same as Modified Euler |
| **RK4** | 0.1353 | **0.0002%** | **Exceptional accuracy** |

**Analytical**: y(2) = e⁻² ≈ 0.1353

## Convergence Analysis

### Theoretical vs Observed Convergence Rates

The convergence analysis perfectly demonstrates theoretical behavior:

| h | Euler Error | RK2 Error | RK4 Error |
|---|-------------|-----------|-----------|
| 0.10 | 7.45e+00 | 9.84e-01 | 4.27e-03 |
| 0.05 | 4.51e+00 | 2.84e-01 | 3.06e-04 |
| 0.025 | 2.52e+00 | 7.63e-02 | 2.05e-05 |
| 0.0125 | 1.34e+00 | 1.97e-02 | 1.33e-06 |

**Error Reduction Factors** (when h halves):
- **Euler**: ~2× (O(h¹) confirmed)
- **RK2**: ~4× (O(h²) confirmed)
- **RK4**: ~16× (O(h⁴) confirmed)

## Algorithm Complexity

| Method | Time Complexity | Function Evaluations | Order of Accuracy |
|--------|----------------|---------------------|-------------------|
| Euler | O(n) | 1 per step | O(h¹) |
| Modified Euler | O(n) | 2 per step | O(h²) |
| RK2 | O(n) | 2 per step | O(h²) |
| RK4 | O(n) | 4 per step | O(h⁴) |

Where n = (x_end - x_start) / h is the number of steps.

## Error Analysis

### Sources of Error
1. **Truncation Error**: From finite difference approximation
2. **Round-off Error**: From floating-point arithmetic
3. **Propagation Error**: Accumulated from previous steps

### Error Estimation
- **Richardson Extrapolation**: Compare solutions with h and h/2
- **Theoretical Bounds**: Based on method order and step size
- **Empirical Validation**: Through analytical comparison

## Applications

### Scientific Computing
- **Population Dynamics**: Growth and decay models
- **Chemical Kinetics**: Reaction rate equations
- **Physics Simulations**: Motion under various forces

### Engineering Applications
- **Control Systems**: System response analysis
- **Circuit Analysis**: RC, RL, RLC circuit behavior
- **Mechanical Systems**: Oscillations and vibrations

### Mathematical Modeling
- **Epidemiology**: Disease spread models (SIR, SEIR)
- **Economics**: Market dynamics and growth models
- **Environmental Science**: Pollution dispersion models

## Method Selection Guidelines

### Choose **Euler's Method** when:
- Simple implementation needed
- Educational or debugging purposes
- Very small step sizes acceptable

### Choose **Modified Euler** when:
- Better accuracy than Euler required
- Computational efficiency important
- Second-order accuracy sufficient

### Choose **RK2** when:
- Good balance of accuracy and speed needed
- Moderate precision requirements
- Memory constraints exist

### Choose **RK4** when:
- High accuracy required
- Standard industrial applications
- Fourth-order precision needed

## Performance Characteristics

### Accuracy Hierarchy (for smooth functions):
1. **RK4**: Superior accuracy (0.0002-0.021% relative error)
2. **Modified Euler/RK2**: Good accuracy (0.36-4.9% relative error)
3. **Euler**: Basic accuracy (10-37% relative error)

### Computational Cost:
- **Function Evaluations per Step**: Euler(1) < Modified Euler/RK2(2) < RK4(4)
- **Memory Usage**: All methods use O(1) additional memory
- **Implementation Complexity**: Euler < Modified Euler ≈ RK2 < RK4

## Stability Analysis

### Stability Regions
- **Euler**: Limited stability region
- **RK2**: Improved stability over Euler
- **RK4**: Excellent stability properties

### Step Size Selection
- **Safety Factor**: Typically 0.8-0.9 of theoretical maximum
- **Adaptive Methods**: Automatic step size adjustment
- **Problem-Dependent**: Stiff equations require smaller steps

## Error Handling

The implementation provides comprehensive error handling:

- **Parameter Validation**: Finite values, positive step sizes
- **Convergence Monitoring**: Maximum iteration limits
- **Numerical Stability**: Overflow and underflow detection
- **User Input Validation**: Interactive mode parameter checking

## Limitations

### Method Limitations
- **Stiff Equations**: Explicit methods may require very small steps
- **Discontinuities**: Reduced accuracy at discontinuous points
- **Long Integration Intervals**: Error accumulation over time

### Implementation Limitations
- **Fixed Step Size**: No automatic adaptation (except adaptive variant)
- **Single Variable**: Currently limited to first-order systems
- **Memory Usage**: Full trajectory storage for complete solutions

## Future Enhancements

- **Adaptive Step Size**: Automatic h adjustment based on error estimates
- **Stiff Equation Solvers**: Implicit methods (BDF, Rosenbrock)
- **Higher-Order Systems**: Extension to systems of ODEs
- **Parallel Processing**: Multi-threaded evaluation for large systems
- **Specialized Methods**: Symplectic integrators for Hamiltonian systems

## Validation and Testing

### Test Problem Verification
All test problems include analytical solutions for accuracy verification:
1. **Separable ODEs**: y' = 2y(x+1)
2. **Linear ODEs**: y' = x + y
3. **Exponential Models**: y' = -y

### Convergence Testing
- Systematic step size reduction
- Error rate confirmation
- Stability boundary determination

