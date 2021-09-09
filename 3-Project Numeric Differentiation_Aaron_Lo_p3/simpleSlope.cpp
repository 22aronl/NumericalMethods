#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

typedef pair<float, float> point;



ofstream file;

int numIndices;
int* indices;
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

void computePoints(float lowerbound, float higherbound, int num)
{
    numIndices = num;
    indices = new int[num];
    float delta = (higherbound - lowerbound) / numIndices;
    for(int i = 0; i < num; i++)
        indices[i] = delta * i + lowerbound;
}

void calculateValuesForArray()
{
    values = new point[numIndices];
    for(int i = 0; i <= numIndices; i++)
        values[i] = point(indices[i], calculateF1(i));
}

int main()
{
    computePoints(-10.0f, 10.0f, 100);
    calculateValuesForArray();
    afdlafj;


    
}