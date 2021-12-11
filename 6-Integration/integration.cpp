/**
 * @file integration.cpp
 * @author Aaron Lo
 * @brief This tests the different integration methods 
 * @version 0.1
 * @date 2021-11-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>
#include <vector>

// define the function to integrate
typedef double (*fn) (double x);
double* xAr;
double** yAr;
double** yArSimpson;
int* histogram;

using namespace std;

/**
 * @brief this defines a function y = x^2 + 1
 * 
 * @param x the x value
 * @return double the y value
 */
double func1(double x) {
    return 1 + x*x;
}

/**
 * @brief this integrates the function y = x^2 + 1 of func 1 into y = x + x^3 / 3
 * 
 * @param x the x value
 * @return double the y value
 */
double integrateFunc1(double x) {
    return x + x * x * x / 3;
}

/**
 * @brief the function y = x * e^-x^2
 * 
 * @param x the x value
 * @return double the y value
 */
double func2(double x) {
    return x * exp(-x*x);
}

/**
 * @brief integrates the function2 y = x * e^-x^2 to y = e^-x^2 / 2
 * 
 * @param x the x value
 * @return double the y value
 */
double integrateFunc2(double x) {
    return -exp(-x*x) / 2;
}

/**
 * @brief function y = x * e^-x
 * 
 * @param x the x value
 * @return double the y value
 */
double func3(double x) {
    return x * exp(-x);
}

/**
 * @brief integrates func3 y = x * e^-x to y = -x * e^-x - e^-x
 * 
 * @param x the x value
 * @return double the y value
 */
double integrateFunc3(double x) {
    return -x*exp(-x) - exp(-x);
}

/**
 * @brief funciton y = sin(x)
 * 
 * @param x the x value
 * @return double the y value
 */
double func4(double x) {
    return sin(x);
}

/**
 * @brief the integrate func4 y = sin(x) to y = - cos(x)
 * 
 * @param x the x value
 * @return double the y value
 */
double integrateFunc4(double x) {
    return -cos(x);
}

/**
 * @brief calculates the integral of a function using the trapezoid rule from the lower bound to the higher bound
 * in iteration number of steps from funcion func
 * 
 * @param lowerBound the lower starting point of hte integral
 * @param higherBound the higher starting point of the integral
 * @param iterations the number of iterations between hte lower and high bound
 * @param func the funcion being integrated
 * @return double the integral of the function
 */
double trapezoidalRule(double lowerBound, double higherBound, double iterations, fn func) {
    double h = (higherBound - lowerBound) / iterations;
    double sum = 0;
    for(int i = 0; i < iterations; i++) {
        sum += h/2 * (func(lowerBound + i*h) + func(lowerBound + (i+1)*h));
    }
    return sum;
}

/**
 * @brief this calculates the integral of a function using the simpson rule from the lower bound to the higher bound
 * in iteration number of steps from funcion func
 * 
 * @param lowerBound the lower starting point of hte integral
 * @param higherBound the higher starting point of the integral
 * @param interations the number of iterations between hte lower and high bound
 * @param func the funcion being integrated
 * @return double the integral of the function
 */
double simpsonRule(double lowerBound, double higherBound, double interations, fn func)
{
    double h = (higherBound - lowerBound) / interations;
    double sum = 0;
    for(int i = 0; i < interations; i++) {
        double l = lowerBound + i*h;
        double r = lowerBound + (i+1)*h;
        sum += (h)/6 * (func(l) + 4*func((l + r)/2) + func(r));
    }
    return sum;
}

/**
 * @brief outputs the results of the integration to a file
 * 
 * @param filename the filename of hte file
 * @param xAr the x coords
 * @param yAr the y coords
 * @param size the size of the x and y arrays
 * @param iterations the number of iterations used in the arrays
 */
void outputToFile(string filename, double* xAr, double** yAr, int size, int iterations) {
    ofstream file;
    file.open(filename + ".txt");
    if (!file.is_open()) {
        cout << "File "  << " not found" << endl;
    }

    for(int i = 0; i < iterations; i++) {
        file << xAr[i] << "\t";
        for(int j = 0; j < size; j++)
        {
            file << yAr[j][i] << "\t";
        }
        file << endl;
    }
    file.close();
}

/**
 * @brief calcualtes hte rms error of the integration at the int func
 * 
 * @param xAr the x coords
 * @param yAr the y coords
 * @param func the function being integrated index of the arrays
 * @param funcToIntegrate the real function
 * @param iterations the iterations 
 * @return double the rms error
 */
double calculateRMS(double* xAr, double** yAr, int func, fn funcToIntegrate, int iterations) {
    double sum = 0;
    for(int i = 0; i < iterations; i++) {
        sum += (yAr[func][i] - (funcToIntegrate(xAr[i]) - funcToIntegrate(xAr[0]))) * (yAr[func][i] - (funcToIntegrate(xAr[i])-funcToIntegrate(xAr[0])));
    }
    return sqrt(sum / iterations);
}

/**
 * @brief charaterizes the histogram of the data
 * 
 * @param size the size of the data, number of bins
 * @param low the low bound of data
 * @param high the high bound of data
 * @param input single input data
 */
void characterizeIntoHistorgram(int size, double input, int low, int high)
{
    double dif = (double)(high - low) / size;
    histogram[(int) ((input - low) / dif)]++;
}

/**
 * @brief calcualtes hte area of a dimension sphere of radius 1
 * 
 * @param length the radius of circle
 * @param iterations the number to tries 
 * @param bins the number of bins
 * @param dimension the dimensions of the sphere (number of diemsions)
 * @return double the area of the sphere
 */
double monteCarloArea(double length, int iterations, int bins, int dimension, string fileN)
{   
    ofstream file;
    file.open(fileN + ".txt");
    if (!file.is_open()) {
        cout << "File "  << " not found" << endl;
    }

    ofstream file2;
    file2.open(fileN + "2.txt");
    if (!file2.is_open()) {
        cout << "File "  << " not found" << endl;
    }
    double inCircle = 0;
    double temp, total;
    int ints = 1000;
    for(int j = 0; j < iterations/ints; j++)
    {

        for(int i = 0; i < ints; i++) {
            total = 0;
            for(int j = 0; j < dimension; j++) {
                temp = (double)rand() / RAND_MAX * (length * 2) - length;
                //characterizeIntoHistorgram(bins, temp + length, 0, 2 * length);
                
                total += temp * temp;
            }

            if(total < length * length) {
                inCircle++;
            }


            
        }
        file << j * ints << "\t" << (inCircle/j/ints) << endl;
        file2 << log(j * ints) << "\t" << log(abs(3.14159265359/6 - (inCircle/j/ints))) << endl;
    }
    file.close();
    file2.close();
    return inCircle / iterations;
}



/**
 * @brief this runs all the code for the project
 * 
 * @return int success!
 */
int main() 
{
    srand(time(NULL));
    double lowerBound = -2.0;
    double higherBound = 2.0;
    int iterations = 100;
    int iterations2 = 100;

    fn funcAr[] = {func1, func2, func3, func4};
    fn funcIntegrateAr[] = {integrateFunc1, integrateFunc2, integrateFunc3, integrateFunc4};
    xAr = new double[iterations2];
    yAr = new double*[4];
    yArSimpson = new double*[4];
    histogram = new int[iterations2];
    for(int i = 0; i < 4; i++)
    {
        yAr[i] = new double[iterations2];
        yArSimpson[i] = new double[iterations2];
    }

    cout << "Trapezoidal: " << trapezoidalRule(lowerBound, higherBound, iterations, func1) << endl;
    cout << "Actual: " << integrateFunc1(higherBound) - integrateFunc1(lowerBound) << endl;

    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < iterations2; j++) {
            xAr[j] = lowerBound + j * (higherBound - lowerBound) / iterations2;
            yAr[i][j] = trapezoidalRule(lowerBound, xAr[j], iterations, funcAr[i]);
            yArSimpson[i][j] = simpsonRule(lowerBound, xAr[j], iterations, funcAr[i]);
        }
    }

    for(int i = 0; i < 4; i++) {
        cout << "TRAP: Function " << i << ": " << calculateRMS(xAr, yAr, i, funcIntegrateAr[i], iterations2) << endl;
    }

    for(int i = 0; i < 4; i++) {
        cout << "\%Difference: " << i << ": " << abs(trapezoidalRule(lowerBound, higherBound, iterations, funcAr[i]) - (funcIntegrateAr[i](higherBound)
            - funcIntegrateAr[i](lowerBound))/(funcIntegrateAr[i](higherBound) - funcIntegrateAr[i](lowerBound)) * 100) << endl;
    }


    for(int i = 0; i < 4; i++) {
        cout << "SIMP: Function " << i << ": " << calculateRMS(xAr, yArSimpson, i, funcIntegrateAr[i], iterations2) << endl;
    }

    for(int i = 0; i < 4; i++) {
        cout << "\%Difference: " << i << ": " << abs(simpsonRule(lowerBound, higherBound, iterations, funcAr[i]) - (funcIntegrateAr[i](higherBound)
            - funcIntegrateAr[i](lowerBound))/(funcIntegrateAr[i](higherBound) - funcIntegrateAr[i](lowerBound)) * 100) << endl;
    }

    
    //cout << "Monte Carlo: " << monteCarloArea(1, 10000000, iterations, 2, "./6-Integration/monteCarlo.txt") << endl;
    cout << "PI / 4: " << M_PI / 4 << endl;

    cout << "Monte Carlo: " << monteCarloArea(1, 1000000, iterations, 3, "./6-Integration/monteCarlo") << endl;
    cout << "PI / 6: " << M_PI / 6 << endl;


    // int iterations3 = 1000;
    // int xArMonte[iterations3];
    // double yArMonte[iterations3];
    // int dist = 100;
    // for(int i = 0; i < iterations3; i++) {
    //     xArMonte[i] = i * dist;
    //     yArMonte[i] = monteCarloArea(dist, xArMonte[i], iterations, 3);
    // }

    outputToFile("./6-Integration/outputTrap", xAr, yAr, 4, iterations2);
    outputToFile("./6-Integration/outputSimp", xAr, yArSimpson, 4, iterations2);

    ofstream file;    
    file.open("./6-Integration/hist.txt");
    if (!file.is_open()) {
        cout << "File "  << " not found" << endl;
    }

    for(int i = 0; i < iterations; i++) {
        file << i << "\t" << histogram[i] << endl;
    }
    file.close();

    
}