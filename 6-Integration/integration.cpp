#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>
#include <vector>

typedef double (*fn) (double x);
double* xAr;
double** yAr;
double** yArSimpson;

using namespace std;


double func1(double x) {
    return 1 + x*x;
}

double integrateFunc1(double x) {
    return x + x * x * x / 3;
}

double func2(double x) {
    return x * exp(-x*x);
}

double integrateFunc2(double x) {
    return -exp(-x*x) / 2;
}

double func3(double x) {
    return x * exp(-x);
}

double integrateFunc3(double x) {
    return -x*exp(-x) - exp(-x);
}

double func4(double x) {
    return sin(x);
}

double integrateFunc4(double x) {
    return -cos(x);
}

double trapezoidalRule(double lowerBound, double higherBound, double iterations, fn func) {
    double h = (higherBound - lowerBound) / iterations;
    double sum = 0;
    for(int i = 0; i < iterations; i++) {
        sum += h/2 * (func(lowerBound + i*h) + func(lowerBound + (i+1)*h));
    }
    return sum;
}

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

double calculateRMS(double* xAr, double** yAr, int func, fn funcToIntegrate, int iterations) {
    double sum = 0;
    for(int i = 0; i < iterations; i++) {
        sum += (yAr[func][i] - (funcToIntegrate(xAr[i]) - funcToIntegrate(xAr[0]))) * (yAr[func][i] - (funcToIntegrate(xAr[i])-funcToIntegrate(xAr[0])));
    }
    return sqrt(sum / iterations);
}

double monteCarloArea(double length, int iterations, int dimension)
{
    double inCircle = 0;
    double temp, total;
    for(int i = 0; i < iterations; i++) {
        total = 0;
        for(int j = 0; j < dimension; j++) {
            temp = (double)rand() / RAND_MAX * (length * 2) - length;
            total += temp * temp;
        }

        if(total <= length * length)
            inCircle++;
    }

    return inCircle / iterations;
}


int main() 
{
    double lowerBound = -2.0;
    double higherBound = 2.0;
    int iterations = 100;
    int iterations2 = 100;

    fn funcAr[] = {func1, func2, func3, func4};
    fn funcIntegrateAr[] = {integrateFunc1, integrateFunc2, integrateFunc3, integrateFunc4};
    xAr = new double[iterations2];
    yAr = new double*[4];
    yArSimpson = new double*[4];

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
            yAr[i][j] = simpsonRule(lowerBound, xAr[j], iterations, funcAr[i]);
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


    cout << "Monte Carlo: " << monteCarloArea(2, 10000000, 2) << endl;
    cout << "PI / 4: " << M_PI / 4 << endl;
    cout << "Monte Carlo: " << monteCarloArea(2, 10000000, 3) << endl;
    cout << "PI / 6: " << M_PI / 6 << endl;

    outputToFile("./6-Integration/outputTrap", xAr, yAr, 4, iterations2);
    outputToFile("./6-Integration/outputSimp", xAr, yArSimpson, 4, iterations2);
}