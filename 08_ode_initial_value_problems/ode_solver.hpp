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
 * @brief Professional ODE solver for initial value problems (Cauchy problems)
 * 
 * This class implements multiple numerical methods for solving first-order
 * ordinary differential equations of the form dy/dx = f(x, y) with initial
 * condition y(x₀) = y₀.
 */

class ODESolver 
{
    public:
        /**
         * @brief Available ODE solving methods
         */
        enum class Method 
        {
            EULER,              // Euler's method (1st order)
            MODIFIED_EULER,     // Modified Euler (Heun's method, 2nd order)
            RUNGE_KUTTA_2,      // 2nd order Runge-Kutta
            RUNGE_KUTTA_4,      // 4th order Runge-Kutta
            ADAPTIVE_RK4        // Adaptive RK4 with error control
        };

        /**
         * @brief Solution point containing x and y values
         */

        struct SolutionPoint 
        {
            double x, y;
            SolutionPoint(double x, double y) : x(x), y(y) {}
        };

        /**
         * @brief Complete solution result with metadata
         */

        struct SolutionResult 
        {
            std::vector<SolutionPoint> points;  // Solution trajectory
            double final_value;                 // y(x_end)
            std::size_t steps_taken;           // Number of integration steps
            double estimated_error;            // Error estimate (for adaptive methods)
            std::string method_name;           // Method used
            bool converged;                    // Convergence status
            
            SolutionResult() : final_value(0.0), steps_taken(0), estimated_error(0.0), 
                            method_name("Unknown"), converged(true) {}
        };

    private:
        static constexpr double DEFAULT_TOLERANCE = 1e-8;
        static constexpr double MIN_STEP_SIZE = 1e-12;
        static constexpr double MAX_STEP_SIZE = 1.0;
        static constexpr std::size_t MAX_STEPS = 1000000;

        /**
         * @brief Validates input parameters
         */

        static void validate_parameters(double x0, double y0, double x_end, double h) 
        {
            if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(x_end) || !std::isfinite(h)) 
            {
                throw std::invalid_argument("All parameters must be finite");
            }
            if (std::abs(x_end - x0) < std::numeric_limits<double>::epsilon()) 
            {
                throw std::invalid_argument("x_end must be different from x0");
            }
            if (h <= 0.0) 
            {
                throw std::invalid_argument("Step size h must be positive");
            }
            if (h > std::abs(x_end - x0)) 
            {
                throw std::invalid_argument("Step size h is too large for the given interval");
            }
        }

    public:
        /**
         * @brief Euler's method for solving dy/dx = f(x, y)
         * @param f The differential equation function f(x, y)
         * @param x0 Initial x value
         * @param y0 Initial y value
         * @param x_end Final x value
         * @param h Step size
         * @return Final y value
         */
        
        [[nodiscard]] static double euler_method(
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            validate_parameters(x0, y0, x_end, h);
            
            const double direction = (x_end > x0) ? 1.0 : -1.0;
            h *= direction;
            
            double x = x0;
            double y = y0;
            
            while ((direction > 0 && x < x_end) || (direction < 0 && x > x_end)) 
            {
                // Adjust step size for final step
                if ((direction > 0 && x + h > x_end) || (direction < 0 && x + h < x_end)) 
                {
                    h = x_end - x;
                }
                
                y += h * f(x, y);
                x += h;
            }
            
            return y;
        }

        /**
         * @brief Modified Euler method (Heun's method)
         */

        [[nodiscard]] static double modified_euler_method(
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            validate_parameters(x0, y0, x_end, h);
            
            const double direction = (x_end > x0) ? 1.0 : -1.0;
            h *= direction;
            
            double x = x0;
            double y = y0;
            
            while ((direction > 0 && x < x_end) || (direction < 0 && x > x_end)) 
            {
                // Adjust step size for final step
                if ((direction > 0 && x + h > x_end) || (direction < 0 && x + h < x_end)) 
                {
                    h = x_end - x;
                }
                
                const double k1 = f(x, y);
                const double k2 = f(x + h, y + h * k1);
                
                y += h * (k1 + k2) / 2.0;
                x += h;
            }
            
            return y;
        }

        /**
         * @brief 2nd order Runge-Kutta method
         */

        [[nodiscard]] static double runge_kutta_2(
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            validate_parameters(x0, y0, x_end, h);
            
            const double direction = (x_end > x0) ? 1.0 : -1.0;
            h *= direction;
            
            double x = x0;
            double y = y0;
            
            while ((direction > 0 && x < x_end) || (direction < 0 && x > x_end)) 
            {
                // Adjust step size for final step
                if ((direction > 0 && x + h > x_end) || (direction < 0 && x + h < x_end)) 
                {
                    h = x_end - x;
                }
                
                const double k1 = f(x, y);
                const double k2 = f(x + h / 2.0, y + h * k1 / 2.0);
                
                y += h * k2;
                x += h;
            }
            
            return y;
        }

        /**
         * @brief 4th order Runge-Kutta method
         */

        [[nodiscard]] static double runge_kutta_4(
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            validate_parameters(x0, y0, x_end, h);
            
            const double direction = (x_end > x0) ? 1.0 : -1.0;
            h *= direction;
            
            double x = x0;
            double y = y0;
            
            while ((direction > 0 && x < x_end) || (direction < 0 && x > x_end)) 
            {
                // Adjust step size for final step
                if ((direction > 0 && x + h > x_end) || (direction < 0 && x + h < x_end)) 
                {
                    h = x_end - x;
                }
                
                const double k1 = f(x, y);
                const double k2 = f(x + h / 2.0, y + h * k1 / 2.0);
                const double k3 = f(x + h / 2.0, y + h * k2 / 2.0);
                const double k4 = f(x + h, y + h * k3);
                
                y += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
                x += h;
            }
            
            return y;
        }

        /**
         * @brief Comprehensive solve method with full trajectory
         */

        [[nodiscard]] static SolutionResult solve_complete(
            Method method,
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            validate_parameters(x0, y0, x_end, h);
            
            SolutionResult result;
            result.method_name = get_method_name(method);
            
            const double direction = (x_end > x0) ? 1.0 : -1.0;
            h *= direction;
            
            double x = x0;
            double y = y0;
            
            // Store initial point
            result.points.emplace_back(x, y);
            
            while ((direction > 0 && x < x_end) || (direction < 0 && x > x_end)) 
            {
                // Adjust step size for final step
                if ((direction > 0 && x + h > x_end) || (direction < 0 && x + h < x_end)) 
                {
                    h = x_end - x;
                }
                
                double y_new;
                switch (method) 
                {
                    case Method::EULER: 
                    {
                        y_new = y + h * f(x, y);
                        break;
                    }
                    case Method::MODIFIED_EULER: 
                    {
                        const double k1 = f(x, y);
                        const double k2 = f(x + h, y + h * k1);
                        y_new = y + h * (k1 + k2) / 2.0;
                        break;
                    }
                    case Method::RUNGE_KUTTA_2: 
                    {
                        const double k1 = f(x, y);
                        const double k2 = f(x + h / 2.0, y + h * k1 / 2.0);
                        y_new = y + h * k2;
                        break;
                    }
                    case Method::RUNGE_KUTTA_4: 
                    {
                        const double k1 = f(x, y);
                        const double k2 = f(x + h / 2.0, y + h * k1 / 2.0);
                        const double k3 = f(x + h / 2.0, y + h * k2 / 2.0);
                        const double k4 = f(x + h, y + h * k3);
                        y_new = y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
                        break;
                    }
                    default:
                        throw std::invalid_argument("Unknown method");
                }
                
                y = y_new;
                x += h;
                result.points.emplace_back(x, y);
                result.steps_taken++;
                
                if (result.steps_taken > MAX_STEPS) 
                {
                    throw std::runtime_error("Maximum number of steps exceeded");
                }
            }
            
            result.final_value = y;
            return result;
        }

        /**
         * @brief Simple solve method returning only final value
         */

        [[nodiscard]] static double solve(
            Method method,
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            switch (method) 
            {
                case Method::EULER:
                    return euler_method(f, x0, y0, x_end, h);
                case Method::MODIFIED_EULER:
                    return modified_euler_method(f, x0, y0, x_end, h);
                case Method::RUNGE_KUTTA_2:
                    return runge_kutta_2(f, x0, y0, x_end, h);
                case Method::RUNGE_KUTTA_4:
                    return runge_kutta_4(f, x0, y0, x_end, h);
                default:
                    throw std::invalid_argument("Unknown method");
            }
        }

        /**
         * @brief Returns method name as string
         */

        [[nodiscard]] static std::string get_method_name(Method method) 
        {
            switch (method) 
            {
                case Method::EULER:
                    return "Euler's Method";
                case Method::MODIFIED_EULER:
                    return "Modified Euler (Heun's Method)";
                case Method::RUNGE_KUTTA_2:
                    return "2nd Order Runge-Kutta";
                case Method::RUNGE_KUTTA_4:
                    return "4th Order Runge-Kutta";
                case Method::ADAPTIVE_RK4:
                    return "Adaptive RK4";
                default:
                    return "Unknown Method";
            }
        }

        /**
         * @brief Returns theoretical order of accuracy
         */

        [[nodiscard]] static std::size_t get_order_of_accuracy(Method method) 
        {
            switch (method) 
            {
                case Method::EULER:
                    return 1;
                case Method::MODIFIED_EULER:
                case Method::RUNGE_KUTTA_2:
                    return 2;
                case Method::RUNGE_KUTTA_4:
                case Method::ADAPTIVE_RK4:
                    return 4;
                default:
                    return 0;
            }
        }

        /**
         * @brief Estimates error using Richardson extrapolation
         */

        [[nodiscard]] static double estimate_error(
            Method method,
            const std::function<double(double, double)>& f,
            double x0, double y0, double x_end, double h
        ) 
        {
            const double y_h = solve(method, f, x0, y0, x_end, h);
            const double y_h2 = solve(method, f, x0, y0, x_end, h / 2.0);
            
            const std::size_t order = get_order_of_accuracy(method);
            const double factor = std::pow(2.0, static_cast<double>(order));
            
            return std::abs(y_h2 - y_h) / (factor - 1.0);
        }
    };

} // namespace numerical_methods
