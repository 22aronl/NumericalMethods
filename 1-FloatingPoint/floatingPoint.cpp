/**
 * @author Aaron Lo
 * @version August 23, 2021
 * 
 * This code runs to see the limits of the compiler and language.
 * This then tests the following:
 * 
 * nonzero / zero
 * zero / zero
 * ±Infinity / zero
 * zero / ±Infinity
 * zero × ±Infinity
 * ±Infinity × ±Infinity
 * ±Infinity / ±Infinity
 * Have your program determine the values of ±MIN and ±MAX by calculating them using the basic two algorithms of:
 * 2 × (MIN / 2) ≠ 0
 * and
 * (MAX × 2) / 2  ≠ Infinity
 * and EPS
 */

#include <iostream>

using namespace std;

/**
 * @brief this runs the main program that has all the code testing floating points
 * 
 * @return int 
 */
int main()
{
    float nonzero = 1.0;
    float zero = 0.0;
    float inf = nonzero / zero;
    float neginf = - 1 * nonzero / zero;
    
    cout << "neginf: " << neginf << endl;
    cout << "nonzero / zero: " << nonzero / zero << endl;
    cout << "zero / zero: " << zero / zero << endl;
    cout << "infinity / zero: " << inf / zero << endl;
    cout << "NegInfinity / zero: " << neginf / zero << endl;
    cout << "zero / infinity: " << zero / inf << endl;
    cout << "zero / Neginfinity: " << zero / neginf << endl;
    cout << "zero * infinity: " << zero * inf << endl;
    cout << "zero * Neginfinity: " << zero * neginf << endl;
    cout << "infinity * infinity: " << inf * inf << endl;
    cout << "neginfinity * infinity: " << neginf * inf << endl;
    cout << "neginfinity * neginfinity: " << neginf * neginf << endl;
    cout << "infinity / infinity: " << inf / inf << endl;
    cout << "infinity / neginfinity: " << inf / neginf << endl;
    cout << "neginfinity / neginfinity: " << neginf / neginf << endl;


    cout << "\nDeterming the values of ±MIN and ±MAX by calculating them..." << endl;

    float min = 1.0f;
    while(2 * (min / 2) != 0)
        min /= 2;

    cout << "+- min: " << min << endl;

    float max = 1.0f;
    while((max * 2) / 2 != inf)
        max *= 2;

    cout << "+-max: " << max << "\n\n";

    cout << "Calculating +-EPS" << endl;

    float eps = 1.0f;
    while(1 + (eps/2) != 1)
        eps /= 2;

    cout << "+-eps: " << eps << endl;

    return (0);
}

//g++ -std=c++11 a.cpp