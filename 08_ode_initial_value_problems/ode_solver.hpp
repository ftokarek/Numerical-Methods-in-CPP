#pragma once

#include <vector>
#include <functional>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <algorithm>

namespace numerical_methods {

class ODESolver 
{
    public:
        enum class Method 
        {
            EULER,
            MODIFIED_EULER,
            RUNGE_KUTTA_2,
            RUNGE_KUTTA_4,
            ADAPTIVE_RK4
        };

        struct SolutionPoint 
        {
            double x, y;
            SolutionPoint(double x, double y) : x(x), y(y) {}
        };

        struct SolutionResult 
        {
            std::vector<SolutionPoint> points;
            double final_value;
            std::size_t steps_taken;
            double estimated_error;
            std::string method_name;
            bool converged;
            
            SolutionResult() : final_value(0.0), steps_taken(0), estimated_error(0.0), 
                            method_name("Unknown"), converged(true) {}
        };

    private:
        static constexpr double DEFAULT_TOLERANCE = 1e-8;
        static constexpr double MIN_STEP_SIZE = 1e-12;
        static constexpr double MAX_STEP_SIZE = 1.0;
        static constexpr std::size_t MAX_STEPS = 1000000;

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
                if ((direction > 0 && x + h > x_end) || (direction < 0 && x + h < x_end)) 
                {
                    h = x_end - x;
                }
                
                y += h * f(x, y);
                x += h;
            }
            
            return y;
        }

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
            
            result.points.emplace_back(x, y);
            
            while ((direction > 0 && x < x_end) || (direction < 0 && x > x_end)) 
            {
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

}
