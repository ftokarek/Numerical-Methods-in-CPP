#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "ode_solver.hpp"

using namespace numerical_methods;

struct ODEProblem 
{
    std::function<double(double, double)> f;
    std::function<double(double)> analytical;
    std::string description;
    std::string equation;
    double x0, y0, x_end;
};

std::vector<ODEProblem> get_test_problems() 
{
    return 
    {
        {
            [](double x, double y) { return 2.0 * y * (x + 1.0); },
            [](double x) { return std::exp(x * x + 2.0 * x); },
            "Separable ODE",
            "y' = 2y(x + 1), y(0) = 1",
            0.0, 1.0, 1.0
        },
        {
            [](double x, double y) { return x + y; },
            [](double x) { return -1.0 - x + 1.1 * std::exp(x); },
            "Linear ODE",
            "y' = x + y, y(0) = 0.1",
            0.0, 0.1, 1.0
        },
        {
            [](double x, double y) { return -y; },
            [](double x) { return std::exp(-x); },
            "Exponential decay",
            "y' = -y, y(0) = 1",
            0.0, 1.0, 2.0
        }
    };
}

void display_solution_comparison(const ODEProblem& problem, double h) 
{
    const int width = 20;
    const int precision = 10;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "ODE SOLUTION COMPARISON\n";
    std::cout << std::string(80, '=') << "\n";
    
    std::cout << "\nProblem: " << problem.description << "\n";
    std::cout << "Equation: " << problem.equation << "\n";
    std::cout << "Interval: [" << problem.x0 << ", " << problem.x_end << "]\n";
    std::cout << "Step size: " << h << "\n";
    
    const double analytical = problem.analytical(problem.x_end);
    std::cout << "\nAnalytical solution: " << std::setw(width) << analytical << "\n";
    
    std::cout << "\n" << std::string(80, '-') << "\n";
    std::cout << std::setw(25) << "Method" 
              << std::setw(width) << "Result" 
              << std::setw(15) << "Error"
              << std::setw(15) << "Relative Error (%)"
              << std::setw(10) << "Steps" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    std::vector<ODESolver::Method> methods = {
        ODESolver::Method::EULER,
        ODESolver::Method::MODIFIED_EULER,
        ODESolver::Method::RUNGE_KUTTA_2,
        ODESolver::Method::RUNGE_KUTTA_4
    };
    
    for (auto method : methods) 
    {
        try 
        {
            auto result = ODESolver::solve_complete(method, problem.f, 
                                                   problem.x0, problem.y0, problem.x_end, h);
            
            const double error = std::abs(result.final_value - analytical);
            const double relative_error = (analytical != 0.0) ? 
                (error / std::abs(analytical)) * 100.0 : 0.0;
            
            std::cout << std::setw(25) << ODESolver::get_method_name(method)
                      << std::setw(width) << result.final_value
                      << std::setw(15) << std::scientific << error
                      << std::setw(15) << std::fixed << relative_error
                      << std::setw(10) << result.steps_taken << "\n";
        }
        catch (const std::exception& e) 
        {
            std::cout << std::setw(25) << ODESolver::get_method_name(method)
                      << std::setw(width) << "Error: " << e.what() << "\n";
        }
    }
    
    std::cout << std::string(80, '=') << "\n";
}

void demonstrate_convergence_analysis(const ODEProblem& problem) 
{
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "CONVERGENCE ANALYSIS: " << problem.description << "\n";
    std::cout << std::string(80, '=') << "\n";
    
    const double analytical = problem.analytical(problem.x_end);
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(10) << "h" 
              << std::setw(20) << "Euler Error"
              << std::setw(20) << "RK2 Error"
              << std::setw(20) << "RK4 Error" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (double h : {0.1, 0.05, 0.025, 0.0125, 0.00625}) 
    {
        std::cout << std::setw(10) << h;
        
        for (auto method : {ODESolver::Method::EULER,
                           ODESolver::Method::RUNGE_KUTTA_2,
                           ODESolver::Method::RUNGE_KUTTA_4}) 
        {
            try 
            {
                const double result = ODESolver::solve(method, problem.f, 
                                                     problem.x0, problem.y0, problem.x_end, h);
                const double error = std::abs(result - analytical);
                std::cout << std::setw(20) << std::scientific << error;
            }
            catch (const std::exception&) 
            {
                std::cout << std::setw(20) << "Error";
            }
        }
        std::cout << "\n";
    }
    
    std::cout << std::string(80, '=') << "\n";
}

void interactive_mode() 
{
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "INTERACTIVE MODE\n";
    std::cout << std::string(60, '=') << "\n";
    
    double x0, y0, x_end, h;
    
    std::cout << "\nEnter initial conditions and parameters:\n";
    std::cout << "Initial x (x0): ";
    std::cin >> x0;
    std::cout << "Initial y (y0): ";
    std::cin >> y0;
    std::cout << "Final x (x_end): ";
    std::cin >> x_end;
    std::cout << "Step size (h): ";
    std::cin >> h;
    
    if (std::abs(x_end - x0) < 1e-10) 
    {
        std::cout << "Warning: x_end is too close to x0. Using default values.\n";
        x0 = 0.0; y0 = 1.0; x_end = 1.0; h = 0.1;
    }
    
    if (h > std::abs(x_end - x0)) 
    {
        std::cout << "Warning: Step size too large. Adjusting to reasonable value.\n";
        h = std::abs(x_end - x0) / 10.0;
    }
    
    if (h <= 0) 
    {
        std::cout << "Warning: Invalid step size. Using h = 0.1.\n";
        h = 0.1;
    }
    
    auto test_problem = get_test_problems()[0];
    
    std::cout << "\nSolving: " << test_problem.equation << "\n";
    std::cout << "With your parameters: x0=" << x0 << ", y0=" << y0 
              << ", x_end=" << x_end << ", h=" << h << "\n";
    
    if (std::abs(x_end - x0) > 2.0 || std::abs(y0) > 10.0) 
    {
        std::cout << "\nWarning: Large interval or initial value may lead to extreme results.\n";
    }
    
    try 
    {
        auto result = ODESolver::solve_complete(ODESolver::Method::RUNGE_KUTTA_4, 
                                               test_problem.f, x0, y0, x_end, h);
        
        std::cout << "\nResult using RK4 method:\n";
        std::cout << "Final value: " << std::fixed << std::setprecision(8) << result.final_value << "\n";
        std::cout << "Steps taken: " << result.steps_taken << "\n";
        
        if (std::abs(result.final_value) > 1e10) 
        {
            std::cout << "\nNote: Result is very large due to exponential growth in the ODE.\n";
            std::cout << "This is mathematically correct but may indicate numerical challenges.\n";
        }
    }
    catch (const std::exception& e) 
    {
        std::cout << "Error: " << e.what() << "\n";
    }
}

int main() 
{
    try 
    {
        std::cout << "Professional ODE Solver Suite for Initial Value Problems\n";
        std::cout << std::string(80, '=') << "\n";
        
        auto test_problems = get_test_problems();
        const double h = 0.1;
        
        for (const auto& problem : test_problems) 
        {
            display_solution_comparison(problem, h);
        }
        
        if (!test_problems.empty()) 
        {
            demonstrate_convergence_analysis(test_problems[0]);
        }
        
        char interactive;
        std::cout << "\nWould you like to test custom parameters? (y/n): ";
        std::cin >> interactive;
        
        if (interactive == 'y' || interactive == 'Y') 
        {
            interactive_mode();
        }
        
        return 0;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        return 1;
    }
}
