/**
 * @file findingZeros.cpp 
 * @author Aaron Lo
 * @brief This purpose of the coding assignemnt is to explore finding roots
 * of a functino through the use of two numerical techniques: Bisection and newtown Raphson.
 * Then, it tests these functions through the use of multiple varied testers to determine drawbacks and dissures that might arise
 * 
 * 
 * 
 * @version 0.1
 * @date 2021-10-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>

// Defines fn as a function that takes a double as parameter and returns a double from a functions
typedef double (*fn) (double a);

using namespace std;

/**
 * @brief Calculates a function x^2
 * 
 * @param a input
 * @return double of the funciton 
 */
double func1(double a)
{
    return a * a;
}

/**
 * @brief calculate the function a^2 + 1
 * 
 * @param a input
 * @return double output 
 */
double func2(double a)
{
    return 1 + a * a;
}

/**
 * @brief calculate sfunction e^( - a^2)
 * 
 * @param a input
 * @return double output 
 */
double func3(double a)
{
    return exp(- a * a);
}

/**
 * @brief calcualte function a * e^(-a^2)
 * 
 * @param a input   
 * @return double output
 */
double func4(double a)
{
    return a * exp(- a * a);
}

/**
 * @brief calculates function x * e^(x^2)
 * 
 * @param a input
 * @return double output 
 */
double func5(double a)
{
    return a * exp(a * a);
}

/**
 * @brief calculates the function x * e^(-x)
 * 
 * @param a input
 * @return double output 
 */
double func6(double a)
{
    return a * exp(- a);
}

/**
 * @brief calcualtes the function a * e^a
 * 
 * @param a input
 * @return double output 
 */
double func7(double a)
{
    return a * exp(a);
}

/**
 * @brief calcualtes the function sqrt ( abs(a))
 * 
 * @param a input
 * @return double output 
 */
double func8(double a)
{
    return sqrt(abs(a));
}

/**
 * @brief function for x^3 + 1
 * 
 * @param x input
 * @return double output
 */
double func9(double x)
{
    return x * x * x ;
}

/**
 * @brief this is the bisection method. It works by binary searching the function down to the zero.
 * However, it requres for the given endpoints to be of opposite signs or it will not work.
 * 
 * @param leftEndpoint the left starting point for the bisection method
 * @param rightEndpoint the right starting point for the method
 * @param tolerance wher ethe method will end (how close does the output need to be to 0)
 * @param maxInterations how many times the method will run before giving up
 * @param func the function (which returns a double) that is to be calcualated
 * @return double the x-position; otherwise, if both endpoints y don't have opposite signs, nan will be returned
 */
double bisection(double leftEndpoint, double rightEndpoint, double tolerance, int maxInterations, double (*func) (double))
{
    if((func(leftEndpoint) < 0) == (func(rightEndpoint) < 0))
        return nan("");

    for(int i = 0; i < maxInterations; i++)
    {
        double mid = (rightEndpoint + leftEndpoint) / 2;

        double l = func(leftEndpoint);
        double r = func(rightEndpoint);
        double m = func(mid);
        //cout << l << " " << m << " " << r << " " << tolerance << endl;
        if(abs(m) < tolerance)
            return mid;
        
        if((l < 0) == (m > 0))
            rightEndpoint = mid;
        else if((r < 0) == (m > 0))
            leftEndpoint = mid;
        else
            throw invalid_argument("Infinite Loop");
        
    } //for(int i = 0; i < maxInterations; i++)

    return nan("Max Iterations Exceeded");

}

/**
 * @brief Calcualtes the derivative of the given function when given
 * through the five point stencil method
 * 
 * @param func the function given to be calculated
 * @param mid the point where the derivative should be calculated
 * @param stepSize how far the function should step in either direction for each step
 * @return double the approximation of the derivative of the functino
 */
double fivePointStencil(double (*func) (double), double mid, double stepSize)
{
    return (-1 * func(mid + 2 * stepSize) + 8 * func(mid + stepSize) - 8 * func(mid - stepSize) + func(mid - 2 * stepSize))/(12 * stepSize);
}

double doubleFivePointStencil(double (*func) (double), double mid, double stepSize)
{
    return (-1 * fivePointStencil(func, mid + 2 * stepSize, stepSize) + 8 * fivePointStencil(func, mid + stepSize, stepSize) - 8 * fivePointStencil(func, mid - stepSize, stepSize) + fivePointStencil(func, mid - 2 * stepSize, stepSize))/(12 * stepSize);
}

/**
 * @brief the newton rpahson method works by implementing x1 = x0 + f(x0) / f'(x0) and multiple iterations
 * the derivative of the function is found through the use of the method fivePointStencil, which employs the five point stencil in 
 * finding a derivative.
 * 
 * @param func the function the newtonRaphson will evaluate
 * @param initial the initial starting parameter
 * @param stepSize the step size in which the five point stencil will use
 * @param maxIterations the max iterations that will be used in the method, exceeding will result in a nan
 * @param tolerance the tolerance to 0 that is needed
 * @return double the location of x where the func hits 0, nan otherwise if max iteration or the 0 doesn't exist
 */
double newtonRaphson(double (*func) (double), double initial, double stepSize, int maxIterations, double tolerance)
{
    double cur = initial;
    for(int i = 0; i < maxIterations; i++)
    {
        double deriv = fivePointStencil(func, cur, stepSize);
        if(deriv == 0)
        {
            cout << "STATIONARY POINT" << endl;
            return nan("");
        }
        cur = cur - func(cur) / deriv;
        if(abs(func(cur)) < tolerance)
            return cur;

        //cout << "STEP " << cur << " " << func(cur) << endl;
    } //for(int i = 0; i < maxIterations; i++)
    return nan("maxIterations Exceeded");
}

double doubleNewtonRaphson(double (*func) (double), double initial, double stepSize, int maxIterations, double tolerance)
{
    double cur = initial;
    for(int i = 0; i < maxIterations; i++)
    {
        double deriv = doubleFivePointStencil(func, cur, stepSize);
        if(deriv == 0)
        {
            cout << "STATIONARY POINT" << endl;
            return nan("");
        }
        cur = cur - fivePointStencil(func, cur, stepSize) / deriv;
        if(abs(fivePointStencil(func, cur, stepSize)) < tolerance)
            return cur;

        //cout << "STEP " << cur << " " << func(cur) << endl;
    } //for(int i = 0; i < maxIterations; i++)
    return nan("maxIterations Exceeded");
}



/**
 * @brief Runs 9 functions to through the testing of the bisection and newtonRaphson. 
 * The bisections tests from a range from -10 to 10 while newton raphson starts at 1.
 * 
 * @return int 0
 */
int main()
{
    fn ar[] = {
        func1, func2, func3, func4, func5, func6, func7, func8, func9
    }; 
    //fn ar[] = {func9};

    for(int i = 0; i < sizeof(ar) / sizeof (ar[0]); i++)
    {
        double bi = bisection(-10, 10, 0.00001, 10000, ar[i]);
        double newt = newtonRaphson(ar[i],0.2, 0.001, 1000000, 0.00001);
        double localMin = doubleNewtonRaphson(ar[i],0.2, 0.001, 1000000, 0.00001);
        cout << "Funciton " << (i + 1) << endl;
        cout << "Bisection: " << bi << "   Newton: " << newt << "    Local Min: " << localMin << endl;
        cout << "Calcuations f" << (i + 1) << "(x)->  B: " << ar[i](bi) << "  N: " << ar[i](newt) << "  L: " << ar[i](localMin) << endl;
        // cout << "BISEcTION " << bi << "      Calculated bi: " << ar[i](bi) << endl;
        // cout << "Newton " << newt << "      Calculated Newton: " << ar[i](newt) << endl;
        // cout << "Local Min: " << localMin << "      Calculated Min: " << ar[i](localMin) << endl;
        cout << endl;
    } //for(int i = 0; i < sizeof(ar) / sizeof (ar[0]); i++)

    return (0);
}