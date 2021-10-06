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
    return x * x * x + 1;
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
 * @return double the x-position
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
        
    }

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

/**
 * @brief the newton rpahs
 * 
 * @param func 
 * @param initial 
 * @param stepSize 
 * @param maxIterations 
 * @param tolerance 
 * @return double 
 */
double newtonRaphson(double (*func) (double), double initial, double stepSize, int maxIterations, double tolerance)
{
    double cur = initial;
    for(int i = 0; i < maxIterations; i++)
    {
        cur = cur - func(cur) / fivePointStencil(func, cur, stepSize);
        if(abs(func(cur)) < tolerance)
            return cur;

        // cout << "STEP " << cur << " " << func(cur) << endl;
    }
    return nan("maxIterations Exceeded");
}


int main()
{
    fn ar[] = {
        func1, func2, func3, func4, func5, func6, func7, func8, func9
    }; 
    //fn ar[] = {func9};

    for(int i = 0; i < sizeof(ar) / sizeof (ar[0]); i++)
    {
        double bi = bisection(-10, 10, 0.00001, 10000, ar[i]);
        double newt = newtonRaphson(ar[i], 1, 0.01, 1000000, 0.000001);
        cout << "BISEcTION " << i << " " << bi << " Calculated bi: " << ar[i](bi) << endl;
        cout << "Newton " << newt << " Calculated Newton: " << ar[i](newt) << endl;
    }
}