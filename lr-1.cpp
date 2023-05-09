#include <vector>
#include <cmath>
#include <iostream>


double f(int a, double x, double n) {
    return (a / (1 + pow(x, n))) - x;
}

double diff_f(int a, double x, double n) {
    return (a * n * pow(x, n-1))/((1+ pow(x, n)*(1+ pow(x, n)))) - 1;
}

double bisection_method(double a, double left, double right, double eps, double n) {
    double x = (left + right) / 2;
    while (abs(f(a, x, n)) >= eps) {
        x = (left + right) / 2;
        if ((f(a, left, n) * f(a, x, n)) < 0) {
            right = x;
        } else {
            left = x;
        }
    }
    return (left + right) / 2;
}

double newtons_method(double a, double left, double right, double eps, int n) {
    double xkMinus1 = (left + right) / 2;
    double xk = xkMinus1 - (f(a, xkMinus1, n) / diff_f(a, xkMinus1, n));
    while (abs(xk - xkMinus1) >= eps) {
        xkMinus1 = xk;
        xk = xkMinus1 - (f(a, xkMinus1, n) / diff_f(a, xkMinus1, n));
    }
    return (xk + xkMinus1) / 2;
}

int main() {

    std::vector<double> segmLength;

    double left = 0;
    double right = 100;
    double eps = 1e-30;
    double a = 1;
    double result;
    double result1;
    double n = 2;

    result = bisection_method(a, left, right, eps, n);
    result1 = newtons_method(a, left, right, eps, n);

    std::cout << result << std::endl;
    std::cout << result1 << std::endl;
    std::cout << segmLength.size();

    return 0;
}
