#pragma once

#include <vector>
#include <functional>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <algorithm>

namespace numerical_methods {

/**
 * @brief Professional numerical integration class implementing multiple quadrature methods
 * 
 * This class provides implementations of various numerical integration techniques
 * including Rectangle (midpoint), Trapezoidal, and Simpson's rules with adaptive
 * error estimation and comprehensive validation.
 */
class NumericalIntegrator {
public:
    /**
     * @brief Integration methods available
     */
    enum class Method {
        RECTANGLE,    // Midpoint rule
        TRAPEZOIDAL,  // Trapezoidal rule
        SIMPSON,      // Simpson's 1/3 rule
        ADAPTIVE_SIMPSON  // Adaptive Simpson with error control
    };

private:
    static constexpr double DEFAULT_TOLERANCE = 1e-8;
    static constexpr std::size_t MAX_RECURSION_DEPTH = 50;
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e6;

    /**
     * @brief Validates integration parameters
     * @param a Lower bound
     * @param b Upper bound
     * @param n Number of intervals
     * @throws std::invalid_argument for invalid parameters
     */
    static void validate_parameters(double a, double b, std::size_t n) {
        if (!std::isfinite(a) || !std::isfinite(b)) {
            throw std::invalid_argument("Integration bounds must be finite");
        }
        if (std::abs(b - a) < EPSILON) {
            throw std::invalid_argument("Integration interval is too small");
        }
        if (n == 0) {
            throw std::invalid_argument("Number of intervals must be positive");
        }
    }

    /**
     * @brief Recursive adaptive Simpson integration
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param tolerance Error tolerance
     * @param depth Current recursion depth
     * @return Integral approximation
     */
    static double adaptive_simpson_recursive(
        const std::function<double(double)>& f,
        double a, double b, double tolerance, std::size_t depth
    ) {
        if (depth > MAX_RECURSION_DEPTH) {
            throw std::runtime_error("Maximum recursion depth exceeded in adaptive Simpson");
        }

        const double c = (a + b) / 2.0;
        const double h = (b - a) / 6.0;
        
        const double fa = f(a);
        const double fb = f(b);
        const double fc = f(c);
        
        // Simpson's rule for the whole interval
        const double S = h * (fa + 4.0 * fc + fb);
        
        // Simpson's rule for left and right halves
        const double d = (a + c) / 2.0;
        const double e = (c + b) / 2.0;
        const double fd = f(d);
        const double fe = f(e);
        
        const double S_left = (h / 2.0) * (fa + 4.0 * fd + fc);
        const double S_right = (h / 2.0) * (fc + 4.0 * fe + fb);
        const double S_split = S_left + S_right;
        
        // Error estimate
        const double error = std::abs(S_split - S) / 15.0;
        
        if (error <= tolerance) {
            return S_split + (S_split - S) / 15.0; // Richardson extrapolation
        } else {
            return adaptive_simpson_recursive(f, a, c, tolerance / 2.0, depth + 1) +
                   adaptive_simpson_recursive(f, c, b, tolerance / 2.0, depth + 1);
        }
    }

public:
    /**
     * @brief Rectangle (midpoint) rule integration
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param n Number of intervals
     * @return Integral approximation
     */
    [[nodiscard]] static double rectangle_rule(
        const std::function<double(double)>& f,
        double a, double b, std::size_t n
    ) {
        validate_parameters(a, b, n);
        
        const double h = (b - a) / static_cast<double>(n);
        double result = 0.0;
        
        for (std::size_t i = 0; i < n; ++i) {
            const double x_mid = a + (static_cast<double>(i) + 0.5) * h;
            result += f(x_mid);
        }
        
        return result * h;
    }

    /**
     * @brief Trapezoidal rule integration
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param n Number of intervals
     * @return Integral approximation
     */
    [[nodiscard]] static double trapezoidal_rule(
        const std::function<double(double)>& f,
        double a, double b, std::size_t n
    ) {
        validate_parameters(a, b, n);
        
        const double h = (b - a) / static_cast<double>(n);
        double result = (f(a) + f(b)) / 2.0;
        
        for (std::size_t i = 1; i < n; ++i) {
            const double x_i = a + static_cast<double>(i) * h;
            result += f(x_i);
        }
        
        return result * h;
    }

    /**
     * @brief Simpson's 1/3 rule integration
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param n Number of intervals (will be made even if odd)
     * @return Integral approximation
     */
    [[nodiscard]] static double simpson_rule(
        const std::function<double(double)>& f,
        double a, double b, std::size_t n
    ) {
        validate_parameters(a, b, n);
        
        // Ensure n is even for Simpson's rule
        if (n % 2 != 0) {
            ++n;
        }
        
        const double h = (b - a) / static_cast<double>(n);
        double result = f(a) + f(b);
        
        for (std::size_t i = 1; i < n; ++i) {
            const double x_i = a + static_cast<double>(i) * h;
            const double multiplier = (i % 2 == 0) ? 2.0 : 4.0;
            result += multiplier * f(x_i);
        }
        
        return result * h / 3.0;
    }

    /**
     * @brief Adaptive Simpson integration with error control
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param tolerance Error tolerance
     * @return Integral approximation
     */
    [[nodiscard]] static double adaptive_simpson(
        const std::function<double(double)>& f,
        double a, double b, double tolerance = DEFAULT_TOLERANCE
    ) {
        if (!std::isfinite(a) || !std::isfinite(b)) {
            throw std::invalid_argument("Integration bounds must be finite");
        }
        if (std::abs(b - a) < EPSILON) {
            throw std::invalid_argument("Integration interval is too small");
        }
        if (tolerance <= 0.0) {
            throw std::invalid_argument("Tolerance must be positive");
        }
        
        return adaptive_simpson_recursive(f, a, b, tolerance, 0);
    }

    /**
     * @brief General integration method dispatcher
     * @param method Integration method to use
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param n Number of intervals (ignored for adaptive methods)
     * @param tolerance Error tolerance (for adaptive methods)
     * @return Integral approximation
     */
    [[nodiscard]] static double integrate(
        Method method,
        const std::function<double(double)>& f,
        double a, double b,
        std::size_t n = 1000,
        double tolerance = DEFAULT_TOLERANCE
    ) {
        switch (method) {
            case Method::RECTANGLE:
                return rectangle_rule(f, a, b, n);
            case Method::TRAPEZOIDAL:
                return trapezoidal_rule(f, a, b, n);
            case Method::SIMPSON:
                return simpson_rule(f, a, b, n);
            case Method::ADAPTIVE_SIMPSON:
                return adaptive_simpson(f, a, b, tolerance);
            default:
                throw std::invalid_argument("Unknown integration method");
        }
    }

    /**
     * @brief Estimates integration error using Richardson extrapolation
     * @param method Integration method
     * @param f Function to integrate
     * @param a Lower bound
     * @param b Upper bound
     * @param n Base number of intervals
     * @return Error estimate
     */
    [[nodiscard]] static double estimate_error(
        Method method,
        const std::function<double(double)>& f,
        double a, double b, std::size_t n
    ) {
        if (method == Method::ADAPTIVE_SIMPSON) {
            throw std::invalid_argument("Error estimation not applicable for adaptive methods");
        }
        
        const double I_n = integrate(method, f, a, b, n);
        const double I_2n = integrate(method, f, a, b, 2 * n);
        
        // Richardson extrapolation error estimate
        double order = 2.0; // Default order
        if (method == Method::SIMPSON) {
            order = 4.0;
        }
        
        return std::abs(I_2n - I_n) / (std::pow(2.0, order) - 1.0);
    }

    /**
     * @brief Returns the theoretical order of accuracy for a method
     * @param method Integration method
     * @return Order of accuracy
     */
    [[nodiscard]] static std::size_t get_order_of_accuracy(Method method) {
        switch (method) {
            case Method::RECTANGLE:
                return 2;
            case Method::TRAPEZOIDAL:
                return 2;
            case Method::SIMPSON:
            case Method::ADAPTIVE_SIMPSON:
                return 4;
            default:
                return 0;
        }
    }

    /**
     * @brief Returns method name as string
     * @param method Integration method
     * @return Method name
     */
    [[nodiscard]] static std::string get_method_name(Method method) {
        switch (method) {
            case Method::RECTANGLE:
                return "Rectangle (Midpoint) Rule";
            case Method::TRAPEZOIDAL:
                return "Trapezoidal Rule";
            case Method::SIMPSON:
                return "Simpson's 1/3 Rule";
            case Method::ADAPTIVE_SIMPSON:
                return "Adaptive Simpson Rule";
            default:
                return "Unknown Method";
        }
    }
};

} // namespace numerical_methods
