#include <iostream>
#include <fstream>
#include <math.h>
#include <tuple>


using namespace std;

typedef pair<float, float> point;



ofstream file;

int numIndices;
float* indices;
point* values;

void openFile(string filename)
{
    file.open(filename);
}

void closeFile()
{
    file.close();
}

float calculateF1(float t)
{
    return exp(- t * t);
}

float calculateF2(float t)
{
    return -2 * t * exp(- t*t);
}

float calculateF3(float t)
{
    return 1.0f * t * t + 6.0f * t + 9.0f;
}

void computePoints(float lowerbound, float higherbound, int num)
{
    numIndices = num;
    indices = new float[num];
    float delta = (higherbound - lowerbound) / numIndices;
    for(int i = 0; i < num; i++)
        indices[i] = delta * i + lowerbound;
}

void calculateValuesForArray(bool one)
{
    values = new point[numIndices];
    for(int i = 0; i <= numIndices; i++)
        if(!one)
            values[i] = point(indices[i], calculateF1(indices[i]));
        else
            values[i] = point(indices[i], calculateF3(indices[i]));
}

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

point* calculateSecondDerivative(point* ar)
{
    point* ret = new point[numIndices - 1];
    for(int i = 0; i < numIndices - 1; i++)
        ret[i+1] = point(indices[i + 1], (ar[i + 1].second - ar[i].second) / (indices[i + 1] - indices[i]));
    return ret;
}

point* calculateMidDerivative(point* ar)
{
    point* ret = new point[numIndices - 1];
    for(int i = 1; i < numIndices - 1; i++)
        ret[i] = point(indices[i], (ar[i + 1].second - ar[i-1].second)/(indices[i+1] - indices[i-1]));
    return ret;
}

void simpleSlope()
{
    openFile("simpleSlope.out");
    pair<point*, pair<point*, point*> > first = calculateFirstDerivative(values);
    // point* secondDeriv = calculateSecondDerivative(values);
    // point* midDeriv = calculateMidDerivative(values);

    point* firstDeriv = first.first;
    point* secondDeriv = first.second.first;
    point* midDeriv = first.second.second;

    // cout << "VALUES" << endl;
    // for(int i = 0; i< numIndices ;i ++)
    // {
    //     cout << i << " " << values[i].first << " " << values[i].second << endl;
    // }
    // cout << endl;

    //calculateValuesForArray(true);

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
    cout << "The best differential method is the mid derivative" << endl;
    closeFile();
}

void threePointDerivative()
{
    point* midDeriv = calculateMidDerivative(values);
    float diffMid = 0.0f;
    for(int i = 1; i < numIndices - 1; i++)
    {
        diffMid += abs(midDeriv[i].second * midDeriv[i].second - calculateF2(midDeriv[i].first) * calculateF2(midDeriv[i].first));
    }
    cout << "DiffMid" << sqrt(diffMid / (numIndices - 2)) << endl;
}

float calculateDerivativeOfParabola(int two, int t)
{
    float t1 = indices[two -1];
    float t2 = indices[two];
    float t3 = indices[two + 1];
    float y1 = values[two - 1].second;
    float y2 = values[two].second;
    float y3 = values[two + 1].second;

    float a = (-t2 *y1 + t3 *y1 + t1* y2 - t3* y2 - t1 *y3 + t2* y3)/((-t1 + t2)* (t2 - t3)* (-t1 + t3));
    float b = (t2 * t2 * y1 - t3 * t3 *y1 - t1 * t1 *y2 + t3*t3 *y2 + t1*t1 *y3 - t2*t2 *y3) / ((t1 - t2) *(t1 - t3) *(t2 - t3));

    return 2.0f * a * t + b;
}

void functionalFit()
{
    float diffParabola = 0.0f;
    for(int i = 1; i < numIndices - 1; i++)
    {
        float deriv = calculateDerivativeOfParabola(i, indices[i]);
        diffParabola += abs(deriv * deriv - calculateF2(indices[i]) * calculateF2(indices[i]));
    }

    cout << "Functional Fit: " << sqrt(diffParabola / (numIndices - 2)) << endl;
}

void fivePointStencil()
{

}

int main()
{
    computePoints(-10.0f, 10.0f, 100);
    calculateValuesForArray(false);
    simpleSlope();
    threePointDerivative();
    functionalFit();
    return (0);
}