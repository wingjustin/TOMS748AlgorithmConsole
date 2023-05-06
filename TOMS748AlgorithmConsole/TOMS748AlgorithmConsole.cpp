// TOMS748AlgorithmConsole.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//

#include <iostream>
#include <sstream>
#include <math.h>
#include "TOMS748Algorithm.h"

#define PRECISION 50

using namespace std;
using namespace MyMath;

double test_func0(double x) {
    return 2 + x * (7 + x * (-14 + x * 5));
}

double test_func1(double x) {
    return 2 + x * (7 + x * (14 + x * 5));
}

double test_func2(double x) {
    return pow(x - 0.3, 5.0);// *(exp(x - 0.3) - 1) / (exp(x - 0.3) - 1);
}

double test_func3(double x) {
    return sin(10 * x * x) * cos(2 * x * x);
}

double test_func4(double x) {
    return x * (exp(x) + 21) / (exp(x) - 1) - 7;
}

double test_func5(double x) {
    return -pow(x, 1.0 / x) + 0.7;
}

double test_func6(double x) {
    return cos(4.7 * x);
}

double test_func7(double x) {
    return sin(cos(exp(x)));
}

int main()
{
    cout << "TOMS 748 Algorithm" << endl;
    cout << "================================================" << endl;

    double root;
    unsigned int max_iter = numeric_limits<unsigned int>::max();
    unsigned int iter = max_iter;

    ostringstream out;
    out.precision(PRECISION);

    cout << "Test 0 : " << endl;
    cout << "f(x) = 2 + 7x - 14x^2 + 5x^3,  locate root in the interval [-1.5, 0]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func0, -1.5, 0, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 1 : " << endl;
    cout << "f(x) = 2 + 7x + 14x^2 + 5x^3,  locate root in the interval [-3, 0]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func1, -3, 0, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 2 : " << endl;
    cout << "f(x) = (x - 0.3)**5,  locate root in the interval [-100, 100]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func2, -100, 100, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 3 : " << endl;
    cout << "f(x) = sin(10 * x ^ 2) * cos(2 * x ^ 2), locate root in the interval[-0.78, -0.2]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func3, -0.78, -0.2, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 4 : " << endl;
    cout << "f(x) = x * (exp(x) + 21) / (exp(x) - 1) - 7,  locate root in the interval [2.1e-16, 6.8]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func4, 2.1e-16, 6.8, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 5 : " << endl;
    cout << "f(x) = -x^(1/x) + 0.7,  locate root in the interval [0, 1.8e7]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func5, 0, 1.8e7, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 6 : " << endl;
    cout << "f(x) = cos(4.7 * x),  locate root in the interval [1.5, 2]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func6, -1.5, 2, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << "Test 7 : " << endl;
    cout << "f(x) = sin(cos(exp(x))),  locate root in the interval [1, 2]" << endl;

    iter = max_iter;
    root = TOMS748Algorithm::FindRoot(test_func7, 1, 2, iter);

    out.str("");
    out.clear();
    out << fixed << root;
    cout << "x = " << move(out).str() << "..." << endl;
    cout << "iter : " << to_string(iter) << endl;
    cout << endl;

    cout << endl << endl;

    system("pause");
}

// 執行程式: Ctrl + F5 或 [偵錯] > [啟動但不偵錯] 功能表
// 偵錯程式: F5 或 [偵錯] > [啟動偵錯] 功能表

// 開始使用的提示: 
//   1. 使用 [方案總管] 視窗，新增/管理檔案
//   2. 使用 [Team Explorer] 視窗，連線到原始檔控制
//   3. 使用 [輸出] 視窗，參閱組建輸出與其他訊息
//   4. 使用 [錯誤清單] 視窗，檢視錯誤
//   5. 前往 [專案] > [新增項目]，建立新的程式碼檔案，或是前往 [專案] > [新增現有項目]，將現有程式碼檔案新增至專案
//   6. 之後要再次開啟此專案時，請前往 [檔案] > [開啟] > [專案]，然後選取 .sln 檔案
