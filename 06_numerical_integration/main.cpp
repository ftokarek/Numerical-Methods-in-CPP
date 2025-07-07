#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "numerical_integrator.hpp"

using namespace numerical_methods;

/**
 * @brief Test function definitions with analytical integrals for verification
 */
struct TestFunction {
    std::function<double(double)> f;
    std::function<double(double, double)> analytical;
    std::string description;
    std::string formula;
};

/**
 * @brief Collection of test functions with known analytical integrals
 */
std::vector<TestFunction> get_test_functions() {
    return {
        {
            [](double x) { return std::sin(x); },
            [](double a, double b) { return -std::cos(b) + std::cos(a); },
            "Sine function",
            "sin(x)"
        },
        {
            [](double x) { return x * x + 2.0 * x + 5.0; },
            [](double a, double b) { 
                auto integral = [](double x) { return x*x*x/3.0 + x*x + 5.0*x; };
                return integral(b) - integral(a);
            },
            "Polynomial function",
            "x² + 2x + 5"
        },
        {
            [](double x) { return std::exp(x); },
            [](double a, double b) { return std::exp(b) - std::exp(a); },
            "Exponential function",
            "exp(x)"
        },
        {
            [](double x) { return 1.0 / (1.0 + x * x); },
            [](double a, double b) { return std::atan(b) - std::atan(a); },
            "Arctangent derivative",
            "1/(1 + x²)"
        }
    };
}

/**
 * @brief Displays comprehensive integration results
 * @param func Test function
 * @param a Lower bound
 * @param b Upper bound
 * @param n Number of intervals
 */
void display_integration_results(const TestFunction& func, double a, double b, std::size_t n) {
    const int width = 20;
    const int precision = 10;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "NUMERICAL INTEGRATION RESULTS\n";
    std::cout << std::string(80, '=') << "\n";
    
    std::cout << "\nFunction: " << func.description << " - f(x) = " << func.formula << "\n";
    std::cout << "Interval: [" << a << ", " << b << "]\n";
    std::cout << "Number of intervals: " << n << "\n";
    
    // Calculate analytical solution
    const double analytical = func.analytical(a, b);
    std::cout << "\nAnalytical solution: " << std::setw(width) << analytical << "\n";
    
    std::cout << "\n" << std::string(80, '-') << "\n";
    std::cout << std::setw(25) << "Method" 
              << std::setw(width) << "Result" 
              << std::setw(15) << "Error"
              << std::setw(15) << "Relative Error (%)" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    // Test different methods
    std::vector<NumericalIntegrator::Method> methods = {
        NumericalIntegrator::Method::RECTANGLE,
        NumericalIntegrator::Method::TRAPEZOIDAL,
        NumericalIntegrator::Method::SIMPSON
    };
    
    for (auto method : methods) {
        try {
            const double result = NumericalIntegrator::integrate(method, func.f, a, b, n);
            const double error = std::abs(result - analytical);
            const double relative_error = (analytical != 0.0) ? (error / std::abs(analytical)) * 100.0 : 0.0;
            
            std::cout << std::setw(25) << NumericalIntegrator::get_method_name(method)
                      << std::setw(width) << result
                      << std::setw(15) << std::scientific << error
                      << std::setw(15) << std::fixed << relative_error << "\n";
        }
        catch (const std::exception& e) {
            std::cout << std::setw(25) << NumericalIntegrator::get_method_name(method)
                      << std::setw(width) << "Error: " << e.what() << "\n";
        }
    }
    
    // Test adaptive Simpson
    try {
        const double adaptive_result = NumericalIntegrator::adaptive_simpson(func.f, a, b, 1e-10);
        const double adaptive_error = std::abs(adaptive_result - analytical);
        const double adaptive_relative_error = (analytical != 0.0) ? 
            (adaptive_error / std::abs(analytical)) * 100.0 : 0.0;
        
        std::cout << std::setw(25) << "Adaptive Simpson"
                  << std::setw(width) << adaptive_result
                  << std::setw(15) << std::scientific << adaptive_error
                  << std::setw(15) << std::fixed << adaptive_relative_error << "\n";
    }
    catch (const std::exception& e) {
        std::cout << std::setw(25) << "Adaptive Simpson"
                  << std::setw(width) << "Error: " << e.what() << "\n";
    }
    
    std::cout << std::string(80, '=') << "\n";
}

/**
 * @brief Demonstrates convergence analysis
 * @param func Test function
 * @param a Lower bound
 * @param b Upper bound
 */
void demonstrate_convergence(const TestFunction& func, double a, double b) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "CONVERGENCE ANALYSIS: " << func.description << "\n";
    std::cout << std::string(80, '=') << "\n";
    
    const double analytical = func.analytical(a, b);
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(10) << "n" 
              << std::setw(20) << "Rectangle Error"
              << std::setw(20) << "Trapezoidal Error"
              << std::setw(20) << "Simpson Error" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (std::size_t n : {10, 20, 50, 100, 200, 500}) {
        std::cout << std::setw(10) << n;
        
        for (auto method : {NumericalIntegrator::Method::RECTANGLE,
                           NumericalIntegrator::Method::TRAPEZOIDAL,
                           NumericalIntegrator::Method::SIMPSON}) {
            try {
                const double result = NumericalIntegrator::integrate(method, func.f, a, b, n);
                const double error = std::abs(result - analytical);
                std::cout << std::setw(20) << std::scientific << error;
            }
            catch (const std::exception&) {
                std::cout << std::setw(20) << "Error";
            }
        }
        std::cout << "\n";
    }
    
    std::cout << std::string(80, '=') << "\n";
}

int main() {
    try {
        std::cout << "Professional Numerical Integration Suite\n";
        std::cout << std::string(80, '=') << "\n";
        
        auto test_functions = get_test_functions();
        
        // Test parameters
        const double a = 0.5;
        const double b = 2.5;
        const std::size_t n = 100;
        
        // Display results for each test function
        for (const auto& func : test_functions) {
            display_integration_results(func, a, b, n);
        }
        
        // Demonstrate convergence for the first function
        if (!test_functions.empty()) {
            demonstrate_convergence(test_functions[0], a, b);
        }
        
        // Interactive mode
        char continue_testing;
        std::cout << "\nWould you like to test custom parameters? (y/n): ";
        std::cin >> continue_testing;
        
        if (continue_testing == 'y' || continue_testing == 'Y') {
            double custom_a, custom_b;
            std::size_t custom_n;
            
            std::cout << "Enter lower bound: ";
            std::cin >> custom_a;
            std::cout << "Enter upper bound: ";
            std::cin >> custom_b;
            std::cout << "Enter number of intervals: ";
            std::cin >> custom_n;
            
            if (!test_functions.empty()) {
                display_integration_results(test_functions[0], custom_a, custom_b, custom_n);
            }
        }
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        return 1;
    }
}