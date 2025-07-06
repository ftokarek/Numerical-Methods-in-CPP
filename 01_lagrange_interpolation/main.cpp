#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

// Struct representing a single interpolation node (x, y)
struct Node 
{
    double x, y;
    Node(double x, double y) : x(x), y(y) {}
};

// Function to compute the Lagrange interpolation polynomial at a given point
double lagrange(const std::vector<Node> &nodes, double point)
{
    if (nodes.empty())
        return 0.0;
    if (nodes.size() == 1)
        return nodes[0].y;

    double result = 0.0;

    for (std::size_t i = 0; i < nodes.size(); ++i) 
    {
        double li = nodes[i].y;
        for (std::size_t j = 0; j < nodes.size(); ++j) 
        {
            if (i != j) 
            {
                li *= (point - nodes[j].x) / (nodes[i].x - nodes[j].x);
            }
        }
        result += li;
    }

    return result;
}

// Function to print interpolation data and result
void printResults(const std::vector<Node> &nodes, double point, double result) 
{
    std::cout << "Number of interpolation nodes: " << nodes.size() << std::endl;
    std::cout << "\nInterpolation nodes:" << std::endl;

    for (const auto &node : nodes) 
    {
        std::cout << std::setw(4) << node.x << " , " << std::setw(4) << node.y << std::endl;
    }

    std::cout << "\nFor point x = " << point << ", interpolated value f(x) = " << result << std::endl;
}

int main() 
{
    // Task 1: Manual input and interpolation
    std::vector<Node> nodes = {{-4, 5}, {-3, 2}, {1, 5}, {2, 2}};

    double point;
    std::cout << "Task 1:\n\nEnter the point at which to evaluate the interpolation: ";
    std::cin >> point;

    double result = lagrange(nodes, point);
    printResults(nodes, point, result);

    // Task 2: Predefined square root-like table interpolation
    std::vector<Node> nodesSqrt = {{27, 3}, {64, 4}, {125, 5}, {216, 6}};
    double resultSqrt = lagrange(nodesSqrt, 50);

    std::cout << "\nTask 2:" << std::endl;
    std::cout << "Interpolated value f(50) = " << resultSqrt << std::endl << std::endl;

    return 0;
}