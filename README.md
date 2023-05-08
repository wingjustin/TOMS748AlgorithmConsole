# TOMS748AlgorithmConsole

 <h3>Algorithm 748 Enclosing Zeros of Continuous Functions</h3>

<p>
不需對連續函數求導的求根算法, 要求初始區間a b 正負相反
<br><br>
主要有4種算法, 二分法, 割線法, 牛頓二次插值 和 逆三次插值
<br><br>
牛頓二次插值 需要3個不同點, 逆三次插值 需要4個不同點
<br><br>
當前算法若果出現估值超出[a,b]區間或除零等壞情況,需要更換算法, 順序為 逆三次插值, 牛頓二次插值, 割線法, 二分法
<br><br>
每次迭代視情況更換算法加速並保證收斂
<br><br>
<div>references:</div>
<ul>
<li>1. boost.org, \boost\math\tools\toms748_solve.hpp</li>
<li>2. document : Alefeld, Potra and Shi: 1995 Algorithm 748 Enclosing Zeros of Continuous Functions</li>
<li>3. https://github.com/jamadagni/toms748/blob/master/toms748.cpp</li>
</ul>
</p>
