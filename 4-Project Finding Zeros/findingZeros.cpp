
#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>

typedef double (*fn) (double a);

using namespace std;


double func1(double a)
{
    return a * a;
}

double func2(double a)
{
    return 1 + a * a;
}

double func3(double a)
{
    return exp(- a * a);
}

double func4(double a)
{
    return a * exp(- a * a);
}

double func5(double a)
{
    return a * exp(a * a);
}

double func6(double a)
{
    return a * exp(- a);
}

double func7(double a)
{
    return a * exp(a);
}

double func8(double a)
{
    return sqrt(abs(a));
}

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

double fivePointStencil(double (*func) (double), double mid, double stepSize)
{
    return (-1 * func(mid + 2 * stepSize) + 8 * func(mid + stepSize) - 8 * func(mid - stepSize) + func(mid - 2 * stepSize))/(12 * stepSize);
}

double newtonRaphson(double (*func) (double), double initial, double stepSize, int maxIterations, double tolerance)
{
    double cur = initial;
    for(int i = 0; i < maxIterations; i++)
    {
        cur = cur - func(cur) / fivePointStencil(func, cur, stepSize);
        if(abs(func(cur)) < tolerance)
            return cur;
    }
    return nan("maxIterations Exceeded");
}


int main()
{
    fn ar[] = {
        func1, func2, func3, func4, func5, func6, func7, func8
    }; 

    for(int i = 0; i < sizeof(ar) / sizeof (ar[0]); i++)
    {
        double bi = bisection(-10, 10, 0.00001, 10000, ar[i]);
        double newt = newtonRaphson(ar[i], 1, 0.001, 1000000, 0.000001);
        cout << "BISEcTION " << i << " " << bi << " Calculated bi: " << ar[i](bi) << endl;
        cout << "Newton " << newt << " Calculated Newton: " << ar[i](newt) << endl;
    }
}