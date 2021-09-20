/**
 * @file simpleSlope.cpp
 * @author Aaron Lo
 * @brief 
 *  This lab seeks to explore how numerical derivatives can be calculated and what drawbacks exist from the different 
 * differences. 
 * The first method explored was the simple slope Derivative. This simply calculates the slope at each point using a left derivative,
 * right derivative, or mid derivative. Then the RMS, SQRT((SUM xf ^2 - xi^2)/N), is calculated of this by the main function calculateF2
 * 
 * The second method explored is the three point derivative. This is calculataed by testing the slope of the points on either side. 
 * Again, the RMS is found.
 * 
 * The functional fit derivative. The parabola for y = a*t^2 + b*t + c is found for the points (y1, t1), (y2, t2), (y3, t3)
 * and the derivative of that parabola, 2a + b, is found. 
 * 
 * Finally, the five point stencil method was used. A derivative was defined as (-f(x + 2h) + 8f(x + h) - 8f(x-h) + f(x-2h)) / (12h).
 * 
 * 
 * @version 0.1
 * @date 2021-09-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>


using namespace std;

//alias for point: a pair of floats: this is a representatino of a point through a pair cfuncitons
typedef pair<float, float> point;

ofstream file;

int numIndices;
float* indices;
point* values;

/**
 * @brief opens a file in the file variable
 * 
 * @param filename the filename to open
 */
void openFile(string filename)
{
    file.open(filename);
}

/**
 * @brief closes the open file
 * 
 */
void closeFile()
{
    file.close();
}

/**
 * @brief calculates the normal functions of e ^ -t^2
 * 
 * @param t the time
 * @return float the function
 */
float calculateF1(float t)
{
    return exp(- t * t);
}

/**
 * @brief calculates the derivatie of calculateF1
 * 
 * 
 * @param t the time
 * @return float the functions -2 * t * e^-t^2
 */
float calculateF2(float t)
{
    return -2 * t * exp(- t*t);
}

/**
 * @brief computes the poitns from the lower bound to higher bound, for num number of points. The
 * generated points are held in the array indices
 * 
 * @param lowerbound the point to start counting from
 * @param higherbound the upper limit to the points
 * @param num the number of points to be generated
 */
void computePoints(float lowerbound, float higherbound, int num)
{
    numIndices = num;
    indices = new float[num];
    float delta = (higherbound - lowerbound) / numIndices;
    for(int i = 0; i < num; i++)
        indices[i] = delta * i + lowerbound;
}

/**
 * @brief this calculates the points generated in the array "indices" with the funciton
 * defined in calculateF1. These are stored as the typedef point array called values
 * 
 */
void calculateValuesForArray()
{
    values = new point[numIndices];
    for(int i = 0; i <= numIndices; i++)
            values[i] = point(indices[i], calculateF1(indices[i]));
}

/**
 * @brief this calculates the derivate of the points in the ar array with simple limits
 * for left, right and mid slope calcualtions
 * 
 * @param ar the array of values to be calcualted
 * @return pair<point*, pair<point*, point*> > this returns for, left, right, and mid repsectively
 */
pair<point*, pair<point*, point*> > calculateFirstDerivative(point* ar)
{
    
    point* first = new point[numIndices - 1];
    point* second = new point[numIndices - 1];
    point* mid = new point[numIndices - 1];

    for(int i = 0; i < numIndices - 1; i++)
    {
        float num = (ar[i + 1].second - ar[i].second) / (indices[i + 1] - indices[i]);
        first[i] = point(indices[i], num);
        second[i] = point(indices[i+1], num);
        mid[i] = point((indices[i] + indices[i+1])/2, num);
    }
    return make_pair(first, make_pair(second, mid));
}

/**
 * @brief this calcualtes the derivative by getting the slope from before and after a certain point
 * 
 * @param ar the values array
 * @return point* the calculated values
 */
point* calculateMidDerivative(point* ar)
{
    point* ret = new point[numIndices - 1];
    for(int i = 1; i < numIndices - 1; i++)
        ret[i] = point(indices[i], (ar[i + 1].second - ar[i-1].second)/(indices[i+1] - indices[i-1]));
    return ret;
}

/**
 * @brief This main function calculates all the RMS for the simple slopes
 * 
 */
void simpleSlope()
{
    openFile("simpleSlope.out");
    pair<point*, pair<point*, point*> > first = calculateFirstDerivative(values);

    point* firstDeriv = first.first;
    point* secondDeriv = first.second.first;
    point* midDeriv = first.second.second;

    double diffOne = 0;
    double diffTwo = 0;
    double diffMid = 0;

    file << indices[0] << " " << firstDeriv[0].second << " " << "NULL SECOND" << " NULL THIRD" << values[0].second << endl;
    for(int i = 1; i < numIndices - 1; i++)
    {
        diffOne += abs(firstDeriv[i].second * firstDeriv[i].second - calculateF2(firstDeriv[i].first)*calculateF2(firstDeriv[i].first));
        diffTwo += abs(secondDeriv[i].second * secondDeriv[i].second - calculateF2(secondDeriv[i].first)*calculateF2(secondDeriv[i].first));
        diffMid += abs(midDeriv[i].second * midDeriv[i].second - calculateF2(midDeriv[i].first)*calculateF2(midDeriv[i].first));
        file << indices[i] << " " << firstDeriv[i].second << " " << secondDeriv[i].second << " " << midDeriv[i].second << " " << values[i].second << endl;
    }

    
    file << indices[numIndices - 1] << " " << "NULL FIRST " << secondDeriv[numIndices - 1].second << " " << midDeriv[numIndices - 1].second << " " << values[numIndices - 1].second << endl;
    
    diffOne /= (numIndices - 2);
    diffTwo /= (numIndices - 2);
    diffMid /= (numIndices - 2);
    cout << "DiffONE " << sqrt(diffOne) << endl;
    cout << "DiffTwo " << sqrt(diffTwo) << endl;
    cout << "DiffMid " << sqrt(diffMid) << endl;
    cout << "The best differential method is the mid derivative\n" << endl;
    closeFile();
}

/**
 * @brief this is the main function that calculates and evaluates the RMS for the 
 * three point derivative
 * 
 */
void threePointDerivative()
{
    point* midDeriv = calculateMidDerivative(values);
    float diffMid = 0.0f;
    for(int i = 1; i < numIndices - 1; i++)
    {
        diffMid += abs(midDeriv[i].second * midDeriv[i].second - calculateF2(midDeriv[i].first) * calculateF2(midDeriv[i].first));
    }
    cout << "Diff Three Point" << sqrt(diffMid / (numIndices - 2)) << endl;
}

/**
 * @brief This is a helper method to find the derivative through the use of a parabola. This parabola is fitted
 * onto three points, and then the derivitive of that is taken.
 * 
 * @param two the center point index
 * @param t the time to be calculated from the derived derivative
 * @return float the derivative specified at time t
 */
float calculateDerivativeOfParabola(int two, float t)
{
    float t1 = indices[two -1];
    float t2 = indices[two];
    float t3 = indices[two + 1];
    float y1 = values[two - 1].second;
    float y2 = values[two].second;
    float y3 = values[two + 1].second;

    // float a = (-t2 *y1 + t3 *y1 + t1* y2 - t3* y2 - t1 *y3 + t2* y3)/((-t1 + t2)* (t2 - t3)* (-t1 + t3));
    // float b = (t2 * t2 * y1 - t3 * t3 *y1 - t1 * t1 *y2 + t3*t3 *y2 + t1*t1 *y3 - t2*t2 *y3) / ((t1 - t2) *(t1 - t3) *(t2 - t3));
    float a = (t3 *(-y1 + y2) + t2* (y1 - y3) + t1 *(-y2 + y3)) / ((t1 - t2) *(t1 - t3)* (t2 - t3));
    float b = (t3 * t3 *  (y1 - y2) + t1*t1* (y2 - y3) + t2*t2* (-y1 + y3))/((t1 - t2) *(t1 - t3)* (t2 - t3));
    return 2.0f * a * t + b;
}

/**
 * @brief This is the method that calculate shte functional fit through the parabola on each point
 * 
 */
void functionalFit()
{
    float diffParabola = 0.0f;
    float deriv = calculateDerivativeOfParabola(1, indices[0]);
    diffParabola += abs(deriv * deriv - calculateF2(indices[0]) * calculateF2(indices[0]));

    for(int i = 1; i < numIndices - 1; i++)
    {
        deriv = calculateDerivativeOfParabola(i, indices[i]);
        diffParabola += abs(deriv * deriv - calculateF2(indices[i]) * calculateF2(indices[i]));
    }

    deriv = calculateDerivativeOfParabola(numIndices - 2, indices[numIndices - 1]);
    diffParabola += abs(deriv * deriv - calculateF2(indices[numIndices - 1]) * calculateF2(indices[numIndices-1]));

    cout << "Functional Fit: " << sqrt(diffParabola / (numIndices)) << endl;
}

/**
 * @brief Tihs calculates a derivative at time t using the five point stencil, requiring two points ahead and two points behind
 * 
 * 
 * @param time the time whree the derivative is calculated
 * @return float the derivative calculated from the stencil 
 */
float calculateFivePointDifferential(int time)
{
    return (- values[time + 2].second + 8 * values[time+1].second - 8 * values[time-1].second + values[time-2].second) / 
        (12 * (indices[0] - indices[1])); //Change the intervale calculation
}

/**
 * @brief This main method calculates the points derivateive with the five point stencil method. Then
 * the rms of thses values are calculated and printed.
 * 
 */
void fivePointStencil()
{
    float diffStencil = 0.0f;
    for(int i = 2; i < numIndices - 2; i++)
    {
        float diff = calculateFivePointDifferential(i);
        diffStencil += abs(diff * diff - calculateF2(indices[i]) * calculateF2(indices[i]));
    }
    cout << "Five Point Stencil: " << sqrt(diffStencil / (numIndices - 4)) << endl;
}

/**
 * @brief the main method that runs all the code. This first computes 100 points to the whole array and calculates
 * the values for that. Then, the simple slope calcualtion is run, RMS calculated. Similar things are done for the
 * three point derivative, functional fit with parabola calculations, and the five point stencil.
 * 
 * @return int success code
 */
int main()
{
    computePoints(-10.0f, 10.0f, 100);
    calculateValuesForArray();
    simpleSlope();
    threePointDerivative();
    functionalFit();
    fivePointStencil();
    return (0);
}