/**
 * @file leastSquare.cpp
 * @author Aaron Lo
 * @brief this runs the least square algorithm to find the best fit line for the given data. It uses
 * gradient descent with moementum to minimize the error.
 * @version 0.1
 * @date 2021-11-10
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>
#include <vector>

// Defines fn as a function that takes a double as parameter and returns a double from a functions
typedef double (*fn) (double x, double* a);
int size;
int arguments;
double* x;
double* y;

using namespace std;

/**
 * @brief This function runs the gaussian formula to calculate the value of the gaussian function. The 
 * formula is: q0^2 * e^(-(x-q1)^2 / (q2^2)) + q3^2
 * 
 * @param x the x value of the gaussian function
 * @param a an array of size 4 that contains the values of the parameters of the gaussian function
 * @return double the y value of the gaussian function
 */
double func1(double x, double* a) {
    return a[0] * a[0] * exp(-(x - a[1]) * (x - a[1]) / (a[2] * a[2])) + a[3] * a[3];
}

/**
 * @brief This function runs a sin function to calculate the value of the sin function. The formula is:
 * q0 * sin(q1 * x + q2) + q3
 * 
 * @param x the x value of the sin function
 * @param a an array of size 4 that contains the values of the parameters of the sin function
 * @return double the y values of the sin function
 */
double func2(double x, double* a) {
    return a[0] * sin(a[1] * x + a[2]) + a[3];
}


/**
 * @brief This reads the data from the file and stores it in the arrays x and y.
 * 
 * @param fileName the file name of the data
 */
void readFile(string fileName)
{
    ifstream file;

    file.open(fileName);
    cout << "Reading File: " << fileName << endl;
    if (!file.is_open()) {
        cout << "File not found" << endl;
        return;
    }

    vector<double> xV;
    vector<double> yV;

    int i = 0;
    double a, b;

    while (!file.eof()) {
        file >> a >> b;
        xV.push_back(a);
        yV.push_back(b);
        i++;
    }
    size = i;

    x = new double[size];
    y = new double[size];

    for(int i = 0; i < size;i ++)
    {
        x[i] = xV[i];
        y[i] = yV[i];
    }

    file.close();
}

/**
 * @brief Prints the values of the parametersi into a file
 * 
 * @param fileName the fila name of the file
 * @param a the parameters of the function
 * @param arguments the number of parameters
 */
void printToOutputFile(string fileName, double* a, int arguments)
{
    ofstream file;

    file.open(fileName);
    cout << "Printing arguments to file: " << fileName << endl;
    if (!file.is_open()) {
        cout << "File not found" << endl;
        return;
    }

    for (int i = 0; i < arguments; i++)
        file << a[i] << endl;

    file.close();
}

/**
 * @brief calculates the partial derivative of the given function with the five point stencil method
 * 
 * @param j the index of the partial derivative
 * @param x the x value of the function
 * @param h the step size of the five point stencil
 * @param a the parameters of the function
 * @param func the function
 * @return double the partial derivative of parameter j of the function
 */
double fivePointPartialDerivative(int j, double x, double h, double* a, fn func)
{
    double ac = a[j];
    a[j] = ac + 2 * h;
    double a1 = func(x, a);
    a[j] = ac + h;
    double a2 = func(x, a);
    a[j] = ac - h;
    double a3 = func(x, a);
    a[j] = ac - 2 * h;
    double a4 = func(x, a);
    a[j] = ac;
    return (-a1 + 8 * a2 - 8 * a3 + a4) / (12 * h);
}

/**
 * @brief calculates the change for the given parameter of hte function using the fivepoint partial derivative method
 * 
 * @param j the index of the parameter
 * @param a the parameters of the functions
 * @param func the function
 * @param h the step size of the five point stencil
 * @return double the change for the parameter j of the function
 */
double calculateFivePointPartialDerivative(int j, double* a, fn func, double h)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
        sum += (y[i] - func(x[i], a)) * fivePointPartialDerivative(j, x[i], h, a, func);
    return sum;
}

/**
 * @brief calcualtes the error of the function using the formale:
 *  E = 1/2 Sigma((y - f(x))^2)
 * 
 * @param a the parameters of the function
 * @param func the function
 * @return double the error of the function
 */
double calculateError(double* a, fn func)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += (y[i] - func(x[i], a)) * (y[i] - func(x[i], a));
    }
    return sum / 2;
}

/**
 * @brief trains the weights of function using the gradient descent method with momentum implemented
 * 
 * @param a the arguments of the function
 * @param func the functions
 * @param threshold the change in error threshold that is needed to stop the training
 * @param learningFactor learning factor of the gradient descent
 * @param maxIterations the maximum number of iterations before the algorithm gives up
 * @param stepSize the step size of the five point stencil that calculate sthe partial derivatives
 * @param momentum momentum of the gradient descent
 * @return double* the updated array of the parameters
 */
double* train(double* a, fn func, double threshold, double learningFactor, int maxIterations, double stepSize, double momentum)
{
    double* deltaq = new double[size];
    double* lastDeltas = new double[size];

    double previousError = calculateError(a, func) + 10000;
    double error = 0;
    int it = 1e5;
    

    for(int i = 0; i < size; i++)
        lastDeltas[i] = 0;

    for(int k = 0; k < maxIterations / it; k++)
    {
        for(int m = 0; m < it; m++)
        {
            
            for(int i = 0; i < arguments; i++)
                deltaq[i] = learningFactor * calculateFivePointPartialDerivative(i, a, func, stepSize);


            for(int i = 0; i < arguments; i++)
            {
                a[i] += deltaq[i] + momentum * lastDeltas[i];
                lastDeltas[i] = deltaq[i];
            }

            error = calculateError(a, func);
            if(abs(error - previousError) < threshold)
            {
                cout << "Error " << error << " " << previousError << " " << error - previousError << endl;
                cout << "Total Iteration: " << m + k * it << endl;
                return a;
            }
            previousError = error;

            if(isnan(error))
                return a;
        } //for(int m = 0; m < it; m++)

        cout << "Iteration " << k * it << endl;
        cout << "Error " << error << "\nLearning Factor " << learningFactor << endl;
        cout << "Parameters ";
        for(int i = 0; i < arguments; i++)
            cout << a[i] << " ";
        cout << "\n\n";
    } //for(int k = 0; k < maxIterations / it; k++)

    cout << "Max Iterations Exceeded: " << learningFactor << endl;
    return a;
}

/**
 * @brief the method taht runs the program: this pulls in the sample data points from a file (default is input.txt but 
 * it can be defined in the command line) and trains the weights of the function using the gradient descent method.
 * After, it outputs the trained weights to a file (default is output.txt but it can be defined in the command line)
 * 
 * @param argc the number of arguments specified in the command line
 * @param argv the arguments specified in the command line
 * @return int 0 if the program runs successfully
 */
int main(int argc, char* argv[])
{
    string basepath = "./5-Least Squares/";
    string defaultFile = "UCrBHe1.txt";
    string outputFile = "output.txt";

    string fileName = basepath + defaultFile;
    // if (argc > 1)
    //     fileName += argv[1];
    // else
    //     fileName += defaultFile;
    
    int curFunc = 1;

    // if(argc > 2)
    //     curFunc = atoi(argv[2]);

    // if (argc > 3)
    //     outputFile = argv[3];

    // cout << "ARGUMENTS" << argc << endl;

    readFile(fileName);


    fn funcAr[] = {func1, func2};
    double a[] = {1, 10, 2, 3};
    double threshold = 1e-5;
    double learningFactor = 1e-6;
    double stepSize = 1e-6;
    double momentum = 0.9;
    int maxIterations = 1e8;
    
    

    arguments = sizeof(a) / sizeof(a[0]);

    cout << "Starting Training -> Learning Factor: " << learningFactor << ", Max Iterations: " << maxIterations << endl;
    cout << "Beginining Error " << calculateError(a, funcAr[curFunc]) << endl << endl;


    double* a1 = train(a, funcAr[curFunc], threshold, learningFactor, maxIterations, stepSize, momentum);
 
    cout << "a0 = " << a1[0] << endl;
    cout << "a1 = " << a1[1] << endl;
    cout << "a2 = " << a1[2] << endl;
    cout << "a3 = " << a1[3] << endl;
    cout << "Error = " << calculateError(a1, funcAr[curFunc]) << endl;

    printToOutputFile(basepath + outputFile, a1, arguments);
    return 0;
}