#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "gauss_legendre_integrator.hpp"

using namespace numerical_methods;

struct TestFunction {
    std::function<double(double)> f;
    std::function<double(double, double)> analytical;
    std::string description;
    std::string formula;
};

std::vector<TestFunction> get_test_functions() {
    return {
        {
            [](double x) { return std::sin(x); },
            [](double a, double b) { return -std::cos(b) + std::cos(a); },
            "Trigonometric function",
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
            "Rational function",
            "1/(1 + x²)"
        }
    };
}

void display_comparison(const TestFunction& func, double a, double b, std::size_t n = 20) {
    const int width = 25;
    const int precision = 12;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(90, '=') << "\n";
    std::cout << "INTEGRATION METHODS COMPARISON\n";
    std::cout << std::string(90, '=') << "\n";
    
    std::cout << "\nFunction: " << func.description << " - f(x) = " << func.formula << "\n";
    std::cout << "Interval: [" << a << ", " << b << "]\n";
    std::cout << "Intervals for classical methods: " << n << "\n";
    
    const double analytical = func.analytical(a, b);
    std::cout << "\nAnalytical solution: " << std::setw(width) << analytical << "\n";
    
    auto results = GaussLegendreIntegrator::benchmark_methods(func.f, a, b, n);
    
    std::cout << "\n" << std::string(90, '-') << "\n";
    std::cout << std::setw(30) << "Method" 
              << std::setw(width) << "Result" 
              << std::setw(15) << "Error"
              << std::setw(15) << "Rel. Error (%)"
              << std::setw(10) << "Func Evals" << "\n";
    std::cout << std::string(90, '-') << "\n";
    
    for (const auto& result : results) {
        const double error = std::abs(result.value - analytical);
        const double relative_error = (analytical != 0.0) ? (error / std::abs(analytical)) * 100.0 : 0.0;
        
        std::cout << std::setw(30) << result.method_name
                  << std::setw(width) << result.value
                  << std::setw(15) << std::scientific << error
                  << std::setw(15) << std::fixed << relative_error
                  << std::setw(10) << result.function_evaluations << "\n";
    }
    
    std::cout << std::string(90, '=') << "\n";
}

void demonstrate_efficiency_analysis() {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "EFFICIENCY ANALYSIS: Function Evaluations vs Accuracy\n";
    std::cout << std::string(80, '=') << "\n";
    
    auto func = [](double x) { return std::exp(x); };
    auto analytical = [](double a, double b) { return std::exp(b) - std::exp(a); };
    
    const double a = 0.0, b = 1.0;
    const double exact = analytical(a, b);
    
    std::cout << "\nFunction: exp(x), Interval: [" << a << ", " << b << "]\n";
    std::cout << "Analytical result: " << std::fixed << std::setprecision(10) << exact << "\n\n";
    
    std::cout << std::setw(25) << "Method" 
              << std::setw(15) << "Func Evals"
              << std::setw(15) << "Error"
              << std::setw(20) << "Efficiency Ratio" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    std::vector<std::pair<GaussLegendreIntegrator::Method, std::size_t>> test_cases = {
        {GaussLegendreIntegrator::Method::RECTANGLE, 100},
        {GaussLegendreIntegrator::Method::TRAPEZOIDAL, 100},
        {GaussLegendreIntegrator::Method::SIMPSON, 100},
        {GaussLegendreIntegrator::Method::GAUSS_LEGENDRE_2, 1},
        {GaussLegendreIntegrator::Method::GAUSS_LEGENDRE_3, 1},
        {GaussLegendreIntegrator::Method::GAUSS_LEGENDRE_4, 1},
        {GaussLegendreIntegrator::Method::GAUSS_LEGENDRE_5, 1},
        {GaussLegendreIntegrator::Method::GAUSS_LEGENDRE_8, 1}
    };
    
    for (const auto& [method, n] : test_cases) {
        auto result = GaussLegendreIntegrator::integrate_advanced(method, func, a, b, n);
        const double error = std::abs(result.value - exact);
        const double efficiency = error > 0 ? 1.0 / (error * result.function_evaluations) : 
                                 std::numeric_limits<double>::infinity();
        
        std::cout << std::setw(25) << result.method_name
                  << std::setw(15) << result.function_evaluations
                  << std::setw(15) << std::scientific << error
                  << std::setw(20) << std::fixed << efficiency << "\n";
    }
    
    std::cout << std::string(80, '=') << "\n";
}

void interactive_mode() {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "INTERACTIVE TESTING MODE\n";
    std::cout << std::string(60, '=') << "\n";
    
    double a, b;
    std::size_t n;
    
    std::cout << "\nEnter integration bounds:\n";
    std::cout << "Lower bound (a): ";
    std::cin >> a;
    std::cout << "Upper bound (b): ";
    std::cin >> b;
    std::cout << "Number of intervals (for classical methods): ";
    std::cin >> n;
    
    auto test_func = TestFunction{
        [](double x) { return std::sin(x); },
        [](double a, double b) { return -std::cos(b) + std::cos(a); },
        "Sine function",
        "sin(x)"
    };
    
    display_comparison(test_func, a, b, n);
}

int main() {
    try {
        std::cout << "Professional Gauss-Legendre Quadrature Integration Suite\n";
        std::cout << std::string(80, '=') << "\n";
        
        auto test_functions = get_test_functions();
        
        std::vector<std::pair<double, double>> intervals = {
            {0.5, 2.5},
            {0.5, 5.0}
        };
        
        for (std::size_t i = 0; i < std::min(test_functions.size(), intervals.size()); ++i) {
            display_comparison(test_functions[i], intervals[i].first, intervals[i].second);
        }
        
        demonstrate_efficiency_analysis();
        
        char interactive;
        std::cout << "\nWould you like to test custom parameters? (y/n): ";
        std::cin >> interactive;
        
        if (interactive == 'y' || interactive == 'Y') {
            interactive_mode();
        }
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        return 1;
    }
}
