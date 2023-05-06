# TOMS748AlgorithmConsole
 How to Find a root of a function using TOMS748 Algorithm

主要有4種算法, 二分法, 割線法, 牛頓二次插值 和 逆三次插值

牛頓二次插值 需要3個不同點, 逆三次插值 需要4個不同點

以上算法若果出現估值超出[a,b]區間或除零等情況,需要更換算法

視情況更換算法加速並保證收斂

references:
* 1. boost.org, \boost\math\tools\toms748_solve.hpp
* 2. document : Alefeld, Potra and Shi: 1995 Algorithm 748 Enclosing Zeros of Continuous Functions
* 3. https://github.com/jamadagni/toms748/blob/master/toms748.cpp
