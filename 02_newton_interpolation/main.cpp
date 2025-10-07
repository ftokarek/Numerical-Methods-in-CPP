#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <filesystem>
#include "newton_interpolator.hpp"

using namespace numerical_methods;

std::vector<Node> read_data_from_file(const std::string& filename) 
{
    if (!std::filesystem::exists(filename)) 
    {
        throw std::runtime_error("File does not exist: " + filename);
    }
    
    std::ifstream file(filename);

    if (!file.is_open()) 
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    int n;

    if (!(file >> n) || n <= 0) 
    {
        throw std::runtime_error("Invalid number of nodes in file (must be positive integer)");
    }

    std::vector<Node> nodes;
    nodes.reserve(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) 
    {
        double x, y;

        if (!(file >> x >> y)) 
        {
            throw std::runtime_error("Invalid data format at line " + std::to_string(i + 2));
        }
        if (!std::isfinite(x) || !std::isfinite(y)) 
        {
            throw std::runtime_error("Non-finite values not allowed at line " + std::to_string(i + 2));
        }
        nodes.emplace_back(x, y);
    }

    return nodes;
}

double get_evaluation_point() 
{
    double input_point;
    std::cout << "Enter the point at which to evaluate the Newton polynomial: ";
    
    if (!(std::cin >> input_point)) 
    {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        throw std::runtime_error("Invalid input for evaluation point");
    }
    
    if (!std::isfinite(input_point)) 
    {
        throw std::runtime_error("Non-finite evaluation point not allowed");
    }
    
    return input_point;
}

void display_results(const NewtonInterpolator& interpolator, double point, double value) 
{
    const int width = 12;
    const int precision = 6;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << "NEWTON INTERPOLATION RESULTS\n";
    std::cout << std::string(50, '=') << "\n\n";
    
    std::cout << "Number of interpolation nodes: " << interpolator.size() << "\n";
    std::cout << "Polynomial degree: " << interpolator.degree() << "\n\n";
    
    std::cout << "Interpolation nodes:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << std::setw(width) << "x" << std::setw(width) << "f(x)" << "\n";
    std::cout << std::string(30, '-') << "\n";
    
    for (const auto& node : interpolator.get_nodes()) 
    {
        std::cout << std::setw(width) << node.x << std::setw(width) << node.y << "\n";
    }
    
    std::cout << "\nEvaluation:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << "x = " << point << "\n";
    std::cout << "P(x) = " << value << "\n\n";

    std::cout << "Newton coefficients (divided differences):\n";
    std::cout << std::string(40, '-') << "\n";

    const auto& coeffs = interpolator.get_coefficients();

    for (std::size_t i = 0; i < coeffs.size(); ++i) 
    {
        std::cout << "b[" << i << "] = " << std::setw(width) << coeffs[i] << "\n";
    }
    
    std::cout << "\n" << std::string(50, '=') << "\n";
}

int main() 
{
    try 
    {
        std::cout << "Newton Polynomial Interpolation\n";
        std::cout << std::string(50, '=') << "\n\n";
        
        auto nodes = read_data_from_file("data.txt");
        
        std::cout << "Successfully loaded " << nodes.size() << " nodes from data.txt\n\n";
        
        NewtonInterpolator interpolator(std::move(nodes));

        double input_point = get_evaluation_point();

        double value = interpolator.evaluate(input_point);
        
        display_results(interpolator, input_point, value);

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";

        return 1;
    }
}
