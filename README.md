# findresonances 

## 原理
返回在用户指定的频率范围内 $[f_{min},f_{max}]$ ，从复数信号的时间轨迹中提取的谐振频率、衰减常数、Q 因子、振幅和相位。findresonances 脚本命令使用一种称为滤波器对角化的谐波反转方法,通过指数衰减的谐波振荡的叠加来近似时间信号，其形式为

$$
s(t) \approx \sum_{k=1}^{N} A_{k}e^{-i(2\pi f_{k}t - \phi_{k})}e^{-\alpha_{k}t}, \text{ for complex signals}
$$

以下是 $N$ 谐振的数量，每个谐振由四个实值参数表征：
* $f_k$：谐振频率。
* $\alpha_k$：衰减常数，其中 $\alpha_k \ge 0$ 。
  或者，衰减由 $Q$ 因子 $Q_k=\omega_k/2\alpha$ 描述，其中 $\omega_k=2\pi f_k$ 是相应的角频率。
* $A_k$：波幅。
* $\phi_k$：相位。
  
 此外，findresonances 返回一个误差估计值，可用于识别命令报告的虚假谐振。这个估计值是相对置信度的衡量标准，即只有通过比较找到的所有共振的估计值才有意义。如果谐振的误差估计值明显大于其余值，则它很有可能是杂散谐振。

 ### 引用
 [1] MANDELSHTAM V A, TAYLOR H S. Harmonic inversion of time signals and its applications[J/OL]. The Journal of Chemical Physics, 1997, 107(17): 6756-6769. DOI:10.1063/1.475324.

## 例子
对以下信号
 ```
    double f1 = 1.765;
	double alpha1 = 0.005;
	double ampl1 = 1.3;
	double phase1 = 0.4;
	double f2 = 2.345;
	double alpha2 = 0.012;
	double ampl2 = 0.45;
	double phase2 = 1.234;

	vec t_long = linspace(0, 20, 201);
	cx_vec signal_long = ampl1 * exp(-IU*(2.0 * PI * f1 * t_long - phase1)) % exp(-alpha1 * t_long) 
		+ ampl2 * exp(-IU*(2 * PI * f2 * t_long - phase2)) % exp(-alpha2 * t_long);
	findresonances(signal_long, t_long, vec{0, 100}, 200, 1e-5, 1000, 1e-10);

 ```
 使用方法
 ```
 findresonances(
	const cx_vec& complex_signal,
	const vec& time_series,
	const vec& freqWin,
	uword j_basis_count=0,
	double removal_criteria=1e-6,
	uword max_iterations=100,
	double error_threshold = 1e-10);
```

### 结果
|频率|衰减常数|Q因子|振幅|相位|误差估计|
|---|---|---|---|---|---|
1.7650e+00  | 5.0000e-03|   1.1090e+03 |  1.3000e+00  | 4.0000e-01  | 4.2243e-16
2.3450e+00   |1.2000e-02  | 6.1392e+02 |  4.5000e-01 |  1.2340e+00 |  1.5187e-16


## 觉得有意思点个赞🤭🥳🤩