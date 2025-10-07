#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "lagrange_interpolator.hpp"

using namespace numerical_methods;

void print_results(const LagrangeInterpolator& interpolator, double point, double result) 
{
    std::cout << "Number of interpolation nodes: " << interpolator.size() << std::endl;
    std::cout << "\nInterpolation nodes:" << std::endl;

    for (const auto& node : interpolator.get_nodes()) 
    {
        std::cout << std::setw(8) << std::fixed << std::setprecision(3) << node.x << " , " << std::setw(8) << node.y << std::endl;
    }

    std::cout << "\nFor point x = " << std::fixed << std::setprecision(6) << point << ", interpolated value f(x) = " << result << std::endl;
}

int main() 
{
    try 
    {
        LagrangeInterpolator interpolator1({{-4, 5}, {-3, 2}, {1, 5}, {2, 2}});

        double point;
        std::cout << "Task 1:\n\nEnter the point at which to evaluate the interpolation: ";
        std::cin >> point;

        double result = interpolator1.evaluate(point);
        print_results(interpolator1, point, result);

        LagrangeInterpolator interpolator2({{27, 3}, {64, 4}, {125, 5}, {216, 6}});
        double result_sqrt = interpolator2.evaluate(50);

        std::cout << "\nTask 2:" << std::endl;
        std::cout << "Interpolated value f(50) = " << std::fixed << std::setprecision(6) << result_sqrt << std::endl;

        return 0;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}