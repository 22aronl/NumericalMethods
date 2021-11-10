#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>

// Defines fn as a function that takes a double as parameter and returns a double from a functions
typedef double (*fn) (double x, double* a);
int size;
int arguments;
double* x;
double* y;

using namespace std;

//Gausian formula
double func1(double x, double* a) {
    //return a[0] * (x - a[1]) * (x - a[1]) + a[2];
    return a[0] * a[0] * exp(-(x - a[1]) * (x - a[1]) / (a[2] * a[2])) + a[3] * a[3];
}

//Sine formula
double func2(double x, double* a) {
    return a[0] * sin(a[1] * x + a[2]) + a[3];
}

//calculate the partial derivative of func1 in respect to a[0]
double deriv10(double x, double* a)
{
    //return (x - a[1]) * (x - a[1]);
    return 2 * a[0] * exp(-(x - a[1]) * (x - a[1]) / (a[2] * a[2]));
}


//calculate the partial derivative of func1 for a[1] into method name deriv12
double deriv11(double x, double* a) {
    //return -2 * (x - a[1]);
    return 2*a[0] * a[0] * exp(-(x - a[1]) * (x - a[1]) / (a[2] * a[2])) * (x - a[1]) / (a[2] * a[2]);
}

double deriv12(double x, double* a) {
    //return 1;
    return 2 * a[0] * a[0] * exp(-(x - a[1]) * (x - a[1]) / (a[2] * a[2])) * (x - a[1]) * (x - a[1]) / (a[2] * a[2] * a[2]);
}

double deriv13(double x, double* a) {
    //return 1;
    return 2 * a[3];
}


void readFile(string fileName)
{
    ifstream file;

    file.open(fileName);
    if (!file.is_open()) {
        cout << "File not found" << endl;
        return;
    }
    file >> size;

    x = new double[size];
    y = new double[size];


    int i = 0;

    while (!file.eof()) {
        file >> x[i] >> y[i];
        i++;
    }
    //size = i-1;
    file.close();
}

void outputFile(string fileName, double* a, int arguments)
{
    ofstream file;

    file.open(fileName);
    if (!file.is_open()) {
        cout << "File not found" << endl;
        return;
    }

    for (int i = 0; i < arguments; i++) {
        file << a[i] << endl;
    }

    file.close();
}

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

double calculateFivePointPartialDerivative(int j, double* a, fn func, double h)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += (y[i] - func(x[i], a)) * fivePointPartialDerivative(j, x[i], h, a, func);
    }
    return sum;
}


double calculatePartialDerivative(int j, double* a, fn* partialAr, fn func)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += (y[i] - func(x[i], a)) * partialAr[j](x[i], a);
        //cout << "x[i] " << x[i] << endl;
        //cout << partialAr[0](1, a) << endl;
        //cout << "sum " << i << " " << j << " " << sum << " " << (y[i] - func(x[i], a)) << " " << partialAr[j](x[i], a) << endl;
    }
    return sum;
}

double calculateError(double* a, fn func)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += (y[i] - func(x[i], a)) * (y[i] - func(x[i], a));
    }
    return sum / 2;
}


double* train(double* a, fn* partialAr, fn func, double threshold, double learningFactor, int maxIterations, double stepSize)
{
    double* deltaq = new double[size];
    double previousError = calculateError(a, func) + 10000;
    double error = 0;
    int it = 1e4 * 5;
    double* lastDeltas = new double[size];
    double MOMENTUM = 0.9;
    for(int i = 0; i < size; i++)
    {
        lastDeltas[i] = 0;
    }

    for(int k = 0; k < maxIterations / it; k++)
    {
        for(int m = 0; m < it; m++)
        {
            
            for(int i = 0; i < arguments; i++)
            {
                deltaq[i] = learningFactor * calculateFivePointPartialDerivative(i, a, func, stepSize);
                //cout << "DELTA Q" << deltaq[i] << endl;
            }

            for(int i = 0; i < arguments; i++)
            {
                a[i] += deltaq[i] + MOMENTUM * lastDeltas[i];
                lastDeltas[i] = deltaq[i];
                //cout << "Arguments " << i << ": " << a[i] << endl;
            }
            error = calculateError(a, func);
            if(abs(error - previousError) < threshold)
            {
                cout << "Total Iteration :" << m + k * it << endl;
                return a;
            }

            // if(error >= previousError)
            //     learningFactor /= 2;
            // else
            //     learningFactor *= 2;
            previousError = error;
            

            //cout << "ERROR " << error << " learning: " << learningFactor << " a[0] " << a[0] << endl;
            if(isnan(error))
                return a;
        }
        cout << "Iteration" << k * it << endl;
        cout << "Error " << error << "\nLearning Factor " << learningFactor << endl;
        cout << "Parameters" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << endl;
        cout << endl;
    }
    cout << "Max Iterations Exceeded: " << learningFactor << endl;
    return a;
}

int main(int argc, char* argv[])
{
    readFile(argc == 1 ? "./5-Least Squares/UCrBHe1.txt" : argv[1]);
    cout << "READING " << size << endl;

    double a[] = {1, 10, 2, 3};
    arguments = sizeof(a) / sizeof(a[0]);
    fn partialAr[] = {deriv10, deriv11, deriv12, deriv13};
    fn funcAr[] = {func1, func2};
    int curFunc = 1;
    cout << "Starting Error " << calculateError(a, funcAr[curFunc]) << endl << endl;
    double threshold = 1e-10;
    double learningFactor = 1e-5;
    int maxIterations = 1e8;
    double stepSize = 1e-6;

    double* a1 = train(a, partialAr, funcAr[curFunc], threshold, learningFactor, maxIterations, stepSize);
 
    cout << "a0 = " << a1[0] << endl;
    cout << "a1 = " << a1[1] << endl;
    cout << "a2 = " << a1[2] << endl;
    cout << "a3 = " << a1[3] << endl;
    cout << "Error = " << calculateError(a1, funcAr[curFunc]) << endl;

    outputFile(argc == 1 ? "./5-Least Squares/output.txt" : argv[2], a1, arguments);
    return 0;
}