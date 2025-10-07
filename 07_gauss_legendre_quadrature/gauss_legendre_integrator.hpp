#pragma once

#include <vector>
#include <functional>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <array>
#include <algorithm>

namespace numerical_methods {

class GaussLegendreIntegrator {
public:
    enum class Method {
        RECTANGLE,
        TRAPEZOIDAL,
        SIMPSON,
        GAUSS_LEGENDRE_2,
        GAUSS_LEGENDRE_3,
        GAUSS_LEGENDRE_4,
        GAUSS_LEGENDRE_5,
        GAUSS_LEGENDRE_8
    };

    struct IntegrationResult {
        double value;
        double estimated_error;
        std::size_t function_evaluations;
        std::string method_name;
        bool high_precision;
        
        IntegrationResult(double val = 0.0) 
            : value(val), estimated_error(0.0), function_evaluations(0),
              method_name("Unknown"), high_precision(false) {}
    };

private:
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e6;

    class FunctionWrapper {
    private:
        std::function<double(double)> f_;
        mutable std::size_t* counter_;

    public:
        FunctionWrapper(const std::function<double(double)>& f, std::size_t* counter)
            : f_(f), counter_(counter) {}

        double operator()(double x) const {
            if (counter_) ++(*counter_);
            return f_(x);
        }
    };

    struct GaussLegendreData {
        std::vector<double> nodes;
        std::vector<double> weights;
        std::size_t order;
    };

    [[nodiscard]] static GaussLegendreData get_gauss_legendre_data(std::size_t n_points) {
        switch (n_points) {
            case 2:
                return {
                    {-0.5773502691896257, 0.5773502691896257},
                    {1.0, 1.0},
                    6
                };
            case 3:
                return {
                    {-0.7745966692414834, 0.0, 0.7745966692414834},
                    {0.5555555555555556, 0.8888888888888888, 0.5555555555555556},
                    8
                };
            case 4:
                return {
                    {-0.8611363115940526, -0.3399810435848563, 
                     0.3399810435848563, 0.8611363115940526},
                    {0.3478548451374538, 0.6521451548625461, 
                     0.6521451548625461, 0.3478548451374538},
                    10
                };
            case 5:
                return {
                    {-0.9061798459386640, -0.5384693101056831, 0.0,
                     0.5384693101056831, 0.9061798459386640},
                    {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                     0.4786286704993665, 0.2369268850561891},
                    12
                };
            case 8:
                return {
                    {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                     -0.1834346424956498, 0.1834346424956498, 0.5255324099163290,
                     0.7966664774136267, 0.9602898564975363},
                    {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                     0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                     0.2223810344533745, 0.1012285362903763},
                    18
                };
            default:
                throw std::invalid_argument("Unsupported number of Gauss-Legendre points: " + std::to_string(n_points));
        }
    }

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

public:
    [[nodiscard]] static double rectangle_rule(
        const FunctionWrapper& f, double a, double b, std::size_t n
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

    [[nodiscard]] static double trapezoidal_rule(
        const FunctionWrapper& f, double a, double b, std::size_t n
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

    [[nodiscard]] static double simpson_rule(
        const FunctionWrapper& f, double a, double b, std::size_t n
    ) {
        validate_parameters(a, b, n);
        
        if (n % 2 != 0) ++n;
        
        const double h = (b - a) / static_cast<double>(n);
        double result = f(a) + f(b);
        
        for (std::size_t i = 1; i < n; ++i) {
            const double x_i = a + static_cast<double>(i) * h;
            const double multiplier = (i % 2 == 0) ? 2.0 : 4.0;
            result += multiplier * f(x_i);
        }
        
        return result * h / 3.0;
    }

    [[nodiscard]] static double gauss_legendre_rule(
        const FunctionWrapper& f, double a, double b, std::size_t n_points
    ) {
        if (!std::isfinite(a) || !std::isfinite(b)) {
            throw std::invalid_argument("Integration bounds must be finite");
        }
        if (std::abs(b - a) < EPSILON) {
            throw std::invalid_argument("Integration interval is too small");
        }

        const auto data = get_gauss_legendre_data(n_points);
        
        const double scale = (b - a) / 2.0;
        const double shift = (a + b) / 2.0;
        
        double result = 0.0;
        for (std::size_t i = 0; i < data.nodes.size(); ++i) {
            const double x = scale * data.nodes[i] + shift;
            result += data.weights[i] * f(x);
        }
        
        return scale * result;
    }

    [[nodiscard]] static IntegrationResult integrate_advanced(
        Method method,
        const std::function<double(double)>& f,
        double a, double b,
        std::size_t n = 1000
    ) {
        std::size_t eval_count = 0;
        FunctionWrapper wrapper(f, &eval_count);
        
        IntegrationResult result;
        result.method_name = get_method_name(method);
        
        try {
            switch (method) {
                case Method::RECTANGLE:
                    result.value = rectangle_rule(wrapper, a, b, n);
                    break;
                case Method::TRAPEZOIDAL:
                    result.value = trapezoidal_rule(wrapper, a, b, n);
                    break;
                case Method::SIMPSON:
                    result.value = simpson_rule(wrapper, a, b, n);
                    break;
                case Method::GAUSS_LEGENDRE_2:
                    result.value = gauss_legendre_rule(wrapper, a, b, 2);
                    result.high_precision = true;
                    break;
                case Method::GAUSS_LEGENDRE_3:
                    result.value = gauss_legendre_rule(wrapper, a, b, 3);
                    result.high_precision = true;
                    break;
                case Method::GAUSS_LEGENDRE_4:
                    result.value = gauss_legendre_rule(wrapper, a, b, 4);
                    result.high_precision = true;
                    break;
                case Method::GAUSS_LEGENDRE_5:
                    result.value = gauss_legendre_rule(wrapper, a, b, 5);
                    result.high_precision = true;
                    break;
                case Method::GAUSS_LEGENDRE_8:
                    result.value = gauss_legendre_rule(wrapper, a, b, 8);
                    result.high_precision = true;
                    break;
                default:
                    throw std::invalid_argument("Unknown integration method");
            }
        } catch (const std::exception&) {
            result.value = std::numeric_limits<double>::quiet_NaN();
            throw;
        }
        
        result.function_evaluations = eval_count;
        
        if (!result.high_precision && n >= 4) {
            try {
                std::size_t dummy_count = 0;
                FunctionWrapper wrapper2(f, &dummy_count);
                const double result_2n = (method == Method::SIMPSON) ?
                    simpson_rule(wrapper2, a, b, 2 * n) :
                    ((method == Method::TRAPEZOIDAL) ?
                        trapezoidal_rule(wrapper2, a, b, 2 * n) :
                        rectangle_rule(wrapper2, a, b, 2 * n));
                
                const double order = (method == Method::SIMPSON) ? 4.0 : 2.0;
                result.estimated_error = std::abs(result_2n - result.value) / (std::pow(2.0, order) - 1.0);
            } catch (...) {
                result.estimated_error = std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            result.estimated_error = std::pow(10.0, -static_cast<double>(get_order_of_accuracy(method)));
        }
        
        return result;
    }

    [[nodiscard]] static double integrate(
        Method method,
        const std::function<double(double)>& f,
        double a, double b,
        std::size_t n = 1000
    ) {
        return integrate_advanced(method, f, a, b, n).value;
    }

    [[nodiscard]] static std::string get_method_name(Method method) {
        switch (method) {
            case Method::RECTANGLE:
                return "Rectangle Rule";
            case Method::TRAPEZOIDAL:
                return "Trapezoidal Rule";
            case Method::SIMPSON:
                return "Simpson's Rule";
            case Method::GAUSS_LEGENDRE_2:
                return "Gauss-Legendre (2 points)";
            case Method::GAUSS_LEGENDRE_3:
                return "Gauss-Legendre (3 points)";
            case Method::GAUSS_LEGENDRE_4:
                return "Gauss-Legendre (4 points)";
            case Method::GAUSS_LEGENDRE_5:
                return "Gauss-Legendre (5 points)";
            case Method::GAUSS_LEGENDRE_8:
                return "Gauss-Legendre (8 points)";
            default:
                return "Unknown Method";
        }
    }

    [[nodiscard]] static std::size_t get_order_of_accuracy(Method method) {
        switch (method) {
            case Method::RECTANGLE:
            case Method::TRAPEZOIDAL:
                return 2;
            case Method::SIMPSON:
                return 4;
            case Method::GAUSS_LEGENDRE_2:
                return 6;
            case Method::GAUSS_LEGENDRE_3:
                return 8;
            case Method::GAUSS_LEGENDRE_4:
                return 10;
            case Method::GAUSS_LEGENDRE_5:
                return 12;
            case Method::GAUSS_LEGENDRE_8:
                return 18;
            default:
                return 0;
        }
    }

    [[nodiscard]] static std::vector<IntegrationResult> benchmark_methods(
        const std::function<double(double)>& f,
        double a, double b,
        std::size_t n = 20
    ) {
        std::vector<Method> methods = {
            Method::RECTANGLE,
            Method::TRAPEZOIDAL,
            Method::SIMPSON,
            Method::GAUSS_LEGENDRE_2,
            Method::GAUSS_LEGENDRE_3,
            Method::GAUSS_LEGENDRE_4,
            Method::GAUSS_LEGENDRE_5,
            Method::GAUSS_LEGENDRE_8
        };
        
        std::vector<IntegrationResult> results;
        results.reserve(methods.size());
        
        for (auto method : methods) {
            try {
                results.push_back(integrate_advanced(method, f, a, b, n));
            } catch (const std::exception&) {
                IntegrationResult error_result;
                error_result.method_name = get_method_name(method);
                error_result.value = std::numeric_limits<double>::quiet_NaN();
                results.push_back(error_result);
            }
        }
        
        return results;
    }
};

}
