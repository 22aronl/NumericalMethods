/**
 * @author Aaron Lo
 * @version August 23, 2021
 * 
 * This code explores the selective side of selective algorithms, as some algorithms only 
 * work when the given input is large or small. through this code project, the roots formula was explored
 * with new ways of calculation that bypassed float rounding errors. Next, we explored how two mathematically
 * equivalent formulas could give different answers at extremes. Finally, we explored how taylor fractions could be
 * used to represent certain forms that would avoid NAN at smaller values and the creation of selective algorithms to deal with
 * that. Here, I decided the the series would have a size of 6, (they vary as the taylor series gets larger) but due
 * to only needing accuracy around 1.0f and smaller, 6 functions. Then, depending if the input is too small,
 * the taylor series or the extact formula can be used.
 * 
 * 
 * 
 */

#include <iostream>
#include <math.h>

#define taylorSeriousSize 6
//smallestTaylorSin(1.0f)

using namespace std;

/**
 * @brief calculates e^x / (e^x - 1)
 * 
 * @param x the e 
 * @return float the result
 */
float function1(float x)
{
    return exp(x)/(exp(x) - 1);
}

/**
 * @brief calculates 1 / (1-e^-x)
 * 
 * @param x the input
 * @return float the result
 */
float function2(float x)
{
    return 1/(1-exp(-x));
}

/**
 * @brief This returns sin(x + sig) - sin(x)
 * 
 * @param x the target value
 * @param sig a small delta
 * @return float the difference between the two sins
 */
float functionSinSmall(float x, float sig)
{
    return sin(x + sig) - sin(x);
}
/**
 * @brief this calculates sin(x+sig) - sin(x) with a better method of combining the sins with double angle
 * using formula sina - sinb = 2cos(a+b/2)sin(a-b/2)
 * 
 * @param x the target value
 * @param sig a small delta
 * @return float the difference between the two sins with better method of calculating
 */
float betterfunctionSinSmall(float x, float sig)
{
    return 2*cos((2*x + sig)/2)*sin(sig/2);
}

/**
 * @brief this calculates sin(x)/x
 * 
 * @param x the degree to be calculated
 * @return float the intended value
 */
float functionSinc(float x)
{
    return sin(x) / x;
}

/**
 * @brief calculates the power of a float
 * 
 * @param x input
 * @param power the power of the input
 * @return float x^power
 */
float powerX(float x, int power)
{
    float y = 1;
    for(int i = 0; i < power; i++)
        y *= x;
    return y;
}

/**
 * @brief Calculates the factorial
 * 
 * @param x the factioral number
 * @return float x!
 */
float factorial(int x)
{
    float y = 1.0f;
    for(int i = 1; i <= x; i++)
        y *= i;
    return y;
}

/**
 * @brief calculates part of a taylor series for sin (x^y / y!)
 * 
 * @param x the variable 
 * @param y the power
 * @return float the part of taylorseries
 */
float partTaylorSin(float x, int y)
{
    return powerX(x, y)/factorial(y);
}

/**
 * @brief calcualtes part of a taylor series for sinc (x^y-1 / y!)
 * 
 * @param x the variable
 * @param y power
 * @return float the part of the taylorseries
 */
float partTaylorSinc(float x, int y)
{
    return powerX(x, y-1)/factorial(y);
}

/**
 * @brief the taylorserious for sinc up to the yth part
 * 
 * @param x the variable
 * @param y the yth part
 * @return float taylor series
 */
float taylorSeriesSinc(float x, int y)
{
    float sum = 0;
    for(int i = y; i > 0; i--)
    {
        if(i % 2 == 0)
            sum -= partTaylorSinc(x, (2*i-1));
        else
            sum += partTaylorSinc(x, (2*i-1));
    }
    return sum;
}

/**
 * @brief The best sinc function that combines both taylor series and normal calculation
 * 
 * @param x the varaible being ccalculated
 * @return float sin(x)/x
 */
float bestSincFunction(float x)
{
    if(x <= 1E-10f)
        return taylorSeriesSinc(x, taylorSeriousSize);
    return functionSinc(x);
}

/**
 * @brief taylor series of sin with x as variable and y as the number of terms
 * 
 * @param x variable
 * @param y terms
 * @return float taylor serious of sin
 */
float taylorSeriesSin(float x, int y)
{
    float sum = 0;
    for(int i = y; i > 0; i--)
    {
        if(i % 2 == 0)
            sum -= partTaylorSin(x, (2*i-1));
        else
            sum += partTaylorSin(x, (2*i-1));
    }
    return sum;
}

/**
 * @brief This calculates the smaller numerator (number of taylor series parts) needed for the taylor
 * series to equal sin
 * 
 * @param x variable
 * @return int min num of terms needed
 */
int smallestTaylorSin(float x)
{
    int numerator = 1;
    while(taylorSeriesSin(x, numerator) != sin(x))
    {
        if(numerator > 20)
            break;
        numerator++;
        //cout << numerator << " " << taylorSeriesSin(x, numerator) << " " << sin(x) << endl;
        //cout << "difference " << taylorSeriesSin(x, numerator) - sin(x) << endl;
    }
    return numerator;
}


/**
 * @brief This is the main method that runs all the code
 * 
 * @return int 
 */
int main()
{
    cout << "Please enter in a, b, and c:" << endl;
    float a, b, c;
    a = 10;
    b = 10;
    c = 10;
    // cin >> a >> b >> c;

    float root1 = (-b + sqrt(b * b - 4*a*c))/(2*a);
    float root2 = (-b - sqrt(b * b - 4*a*c))/(2*a);

    cout << "x of ax^2+bx+c is: " << root1 << " and " << root2 << endl;

    float x2 = c / (a * root2);

    cout << "x2 is :" << x2 << endl;

    const float large_value = 100000;
    cout << "Evaluating large value " << large_value << " in funciton 1 and 2" << endl;

    cout << "Positive evaluation" << endl;

    cout << "Funciton 1: " << function1(large_value) << " and function 2: " << function2(large_value) << endl;

    cout << "Negative evaluation" << endl;
    cout << "Funciton 1: " << function1(-large_value) << " and function 2: " << function2(-large_value) << endl;

    const float value = 100;
    float delta [] = {1.0f, 0.1f, 0.01f, 0.001f, 0.0001f, 1E-5f, 1E-6f, 1E-7f, 1E-8f};

    for(int i = 0; i < sizeof(delta)/sizeof(delta[0]); i++)
    {
        cout << "Delta Value = " << delta[i] << endl;
        cout << "Normal Value = " << value << endl;
        cout << "Calculated sin(x + delta) - sin(x): " << functionSinSmall(value, delta[i]) << endl;
        cout << "Calculate previous value with a better calculation: " << betterfunctionSinSmall(value, delta[i]) << endl;
        cout << endl;
    } // for(int i = 0; i < sizeof(delta)/sizeof(delta[0]); i++)
    cout << endl;

    float smallValue [] = {5.0f, 1.0f, 0.5f, 0.1f, 1E-2f, 1E-3f, 1E-4f, 1E-5f, 1E-6f, 1E-7f, 1E-8f, 1E-9f, 1E-10f, 0};

    for(int i = 0; i < sizeof(smallValue)/sizeof(smallValue[0]); i++)
    {
        cout << "Small Value: " << smallValue[i] << endl;
        cout << "Normal sin(x)/x calculation: " << functionSinc(smallValue[i]) << endl;
        cout << "Taylor Series sin(x)/x calculation: " << taylorSeriesSinc(smallValue[i], taylorSeriousSize) << endl;
        cout << "Difference: " << taylorSeriesSinc(smallValue[i], taylorSeriousSize) - functionSinc(smallValue[i])<<endl;
        cout << "Best Sinc Calculation: " << bestSincFunction(smallValue[i]) << endl;
        cout << endl;
    } //for(int i = 0; i < sizeof(smallValue)/sizeof(smallValue[0]); i++)

    //cout << "Smallest taylor sin: " << smallestTaylorSin(0.001f) << endl;
    //cout << "SIN: " << sin(1) << endl;
    //cout << taylorSeriesSin(1.0f, 2) << endl;

    return (0);
}