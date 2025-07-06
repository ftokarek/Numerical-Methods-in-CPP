#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

using namespace std;

// Struct representing a single interpolation node (x, y)
struct Node 
{
    double x, y;
    Node(double x, double y) : x(x), y(y) {}
};

// Function to compute the Lagrange interpolation polynomial at a given point
double lagrange(const vector<Node> &nodes, double point)
{
    if (nodes.empty())
        return 0.0;
    if (nodes.size() == 1)
        return nodes[0].y;

    double result = 0.0;

    for (size_t i = 0; i < nodes.size(); ++i) 
    {
        double li = nodes[i].y;
        for (size_t j = 0; j < nodes.size(); ++j) 
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
void printResults(const vector<Node> &nodes, double point, double result) 
{
    cout << "Number of interpolation nodes: " << nodes.size() << endl;
    cout << "\nInterpolation nodes:" << endl;

    for (const auto &node : nodes) 
    {
        cout << setw(4) << node.x << " , " << setw(4) << node.y << endl;
    }

    cout << "\nFor point x = " << point << ", interpolated value f(x) = " << result << endl;
}

int main() 
{
    // Task 1: Manual input and interpolation
    vector<Node> nodes = {{-4, 5}, {-3, 2}, {1, 5}, {2, 2}};

    double point;
    cout << "Task 1:\n\nEnter the point at which to evaluate the interpolation: ";
    cin >> point;

    double result = lagrange(nodes, point);
    printResults(nodes, point, result);

    // Task 2: Predefined square root-like table interpolation
    vector<Node> nodesSqrt = {{27, 3}, {64, 4}, {125, 5}, {216, 6}};
    double resultSqrt = lagrange(nodesSqrt, 50);

    cout << "\nTask 2:" << endl;
    cout << "Interpolated value f(50) = " << resultSqrt << endl << endl;

    return 0;
}