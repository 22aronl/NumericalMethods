#include <iostream>
#include <fstream>
#include <math.h>

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
            values[i] = point(indices[i], calculateF2(indices[i]));
}

point* calculateFirstDerivative(point* ar)
{
    point* ret = new point[numIndices - 1];
    for(int i = 0; i < numIndices - 1; i++)
        ret[i] = point(indices[i], (ar[i + 1].second - ar[i].second) / (indices[i + 1] - indices[i]));
    return ret;
}

point* calculateSecondDerivative(point* ar)
{
    point* ret = new point[numIndices - 1];
    for(int i = 0; i < numIndices - 1; i++)
        ret[i] = point(indices[i + 1], (ar[i + 1].second - ar[i].second) / (indices[i + 1] - indices[i]));
    return ret;
}

point* calculateMidDerivative(point* ar)
{
    point* ret = new point[numIndices - 1];
    for(int i = 1; i < numIndices - 1; i++)
        ret[i] = point(indices[i], (ar[i + 1].second - ar[i-1].second)/(indices[i+1] - indices[i-1]));
    return ret;
}

int main()
{
    openFile("driv.out");
    computePoints(-10.0f, 10.0f, 100);
    calculateValuesForArray(false);
    point* firstDeriv = calculateFirstDerivative(values);
    point* secondDeriv = calculateSecondDerivative(values);
    point* midDeriv = calculateMidDerivative(values);

    // cout << "VALUES" << endl;
    // for(int i = 0; i< numIndices ;i ++)
    // {
    //     cout << i << " " << values[i].first << " " << values[i].second << endl;
    // }
    // cout << endl;

    calculateValuesForArray(true);

    double diffOne = 0;
    double diffTwo = 0;
    double diffMid = 0;

    file << indices[0] << " " << firstDeriv[0].second << " " << "NULL SECOND" << " NULL THIRD" << values[0].second << endl;
    for(int i = 1; i < numIndices - 1; i++)
    {
        diffOne += abs(firstDeriv[i].second * firstDeriv[i].second - values[i].second * values[i].second);
        diffTwo += abs(secondDeriv[i].second * secondDeriv[i].second - values[i].second * values[i].second);
        diffMid += abs(midDeriv[i].second * midDeriv[i].second - values[i].second * values[i].second);
        file << indices[i] << " " << firstDeriv[i].second << " " << secondDeriv[i].second << " " << midDeriv[i].second << " " << values[i].second << endl;
    }

    
    file << indices[numIndices - 1] << " " << "NULL FIRST " << secondDeriv[numIndices - 1].second << " " << midDeriv[numIndices - 1].second << " " << values[numIndices - 1].second << endl;
    
    diffOne /= numIndices;
    diffTwo /= numIndices;
    diffMid /= numIndices;
    cout << "DiffONE " << sqrt(diffOne) << endl;
    cout << "DiffTwo " << sqrt(diffTwo) << endl;
    cout << "DiffMid " << sqrt(diffMid) << endl;
    cout << "The best differential method is the mid derivative" << endl;
    closeFile();
}