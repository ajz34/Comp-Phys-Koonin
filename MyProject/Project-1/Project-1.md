# Project I

## Fundamental Calculation Conceptions

The calculation method is:
\[
\Theta(b) = \int_b^{r_\mathrm{max}} \frac{2b \, \mathrm d r}{r^2} \biggl( 1 - \frac{b^2}{r^2} \biggr)^{-1/2} -
\int_{r_\mathrm{min}}^{r_\mathrm{max}} \frac{2b \, \mathrm d r}{r^2} \biggl( 1 - \frac{b^2}{r^2} - \frac{v(r)}{e} \biggr)^{-1/2}
\]
where
* $r$ and $b$ are the scaled (by Bohr radius $a$) distance and impact factor;
* $e$ is the scaled energy (by the lowest point of the potential $V_0$);
* $v(r)$ is the scaled potential (also by $V_0$); the formula of the Lennard-Jones potential is
\[
v(r) = 4 \big( r^{-12} - r^{-6} \big)
\]
* $r_\mathrm{max}$ is defined manually (scaled by Bohr radius $a$); it should be infinity in principle, however, that is out of our capability to manipulate an integral with it's upper limit be the infinity. It is said to be almost save to set $r_\mathrm{max}$ be $3$; however, I believe it is inadequate in most cases since $b$ can be larger than $r_\mathrm{max}$.
* $r_\mathrm{min}$ is defined by the domain of the following function (which domain is almost the same as the integral core):
\[
f_c(r) = 1 - \frac{b^2}{r^2} - \frac{v(r)}{e}
\]
\[
f_c(r_\mathrm{min}) = 0
\]
* It's very likely to have many roots for the equation, what we want is the greatest root of the equation.

For the FORTRAN code, $e$ is defined at the common block, and won't be called as a variable.

## Mathematica Code for the Fundamental Calculation

This following code shows how the basic calculation of the scattering angle can be calculated:
```Mathematica
e = 0.1; b = 2.258; rmax = 100000;
FuncC[r_] := b*r^-2*(1 - b^2*r^-2 - 4/e*(r^-12 - r^-6))^(-1/2)
```
```Mathematica
rmin = x /.
  NSolve[FuncC[x]^(-1) == 0 && Re[x] > 0 && Abs[Im[x]] < 10^-8, x][[1]]
```
```Mathematica
(2 NIntegrate[b/r^2*(1 - b^2/r^2)^(-1/2), {r, b, rmax}])/
   Pi*180 - (NIntegrate[FuncC[r], {r, rmin, rmax}])/Pi*180
```

## If the Potential be the Square Potential?

To start any discussion, we should first observe the behavior of the integral core $f_c(r)$ as defined before.

These following figures are taken under the situation that $r_\mathrm{max} = 2.0$ (here the $r_\mathrm{min}$ does not means the integral upper limit, but the maximum distance the potential $v(r) = u_0$ can affect).

![](./Figure/Square-Potential-1.png)

The first figure is obtained at $e=-0.5$, and the three curves are obtained at $b= 1.6,2.0,2.4$ (from top to bottom). $x$ axis refers to $r$, while $y$ axis refers to $f_c(r)$. This figure somehow looks like the situation when Lennard-Jones potential is applied and $e$ is small. For $b=2.0$, the electron is "trapped" by the particle and electron orbital merges, as well as $\Theta$ is extremely negative (actually it's a singularity). When $b$ becomes larger or smaller than $2.0$, the value of $\Theta$ becomes larger, and finally converges to $Pi$ (for $b \rightarrow 0$) or $0$ (for $b \rightarrow +\infty$).

![](./Figure/Square-Potential-2.png)

The first figure is obtained at $e=2.0$, and the three curves are obtained at $b= 1.0,1.5,2.0$ (from top to bottom).

Since there are only one possible root for $f_c(r)=0$, the situation is much more simpler than the former situation.

## "Trap" by Central Particle

For a low energy $e$, it is possible to let the electron be "trapped", forming an orbital around the central particle. This situation can happen when there are one minimum for $f_c(r)$ which is slightly **above** zero (but not slightly below zero!). Here's a figure that should explain what happens:

![](./Figure/Trap.png)

This figure is obtained at $e=0.5$, $b=1.91$, Lennard-Jones potential. From this figure, we know that $r_\mathrm{min} \sim 1.19$, so $r=1.54$ should be included by the integral; however, $f_c(1.54)$ is just slightly above zero, so the integral core
\[
\frac{2b}{r^2} \biggl( 1 - \frac{b^2}{r^2} - \frac{v(r)}{e} \biggr)^{-1/2} \rightarrow +\infty
\]
Since the integral of the integral core mentioned above is subtracted by
\[
\pi \sim \int_b^{r_\mathrm{max}} \frac{2b \, \mathrm d r}{r^2} \biggl( 1 - \frac{b^2}{r^2} \biggr)^{-1/2}
\]
so the value of $\Theta$ can be extremely small. The extremely small $\Theta$ implies that the electron have to rotate around the central particle for several times before it leaves the particle and depart to infinity.


## Where is the singularity of the "trap"?

Our problem is that provided a small $e$, find the $b$ as well as $r$ which satisfies the singularity condition; this condition is
\[
\begin{align}
f_c(r) &= 1 - \frac{b^2}{r^2} + \frac 1 e (- 4 r^{-12} + 4 r^{-6}) = 0 \\
f'_ c(r) &= \frac{2 b^2}{r^3} + \frac 1 e (48 r^{-13} - 24 r^{-7}) = 0 \\
f''_ c(r) &= - \frac{b^2}{r^4} + \frac 1 e (-624 r^{-14} + 168 r^{-8}) > 0
\end{align}
\]

First we solve $r$:
\[
f_c(r) + \frac r 2 f'_ c(r) = \frac{e r^{12} - 8 r^6 + 20}{e \, r^{12}} = 0
\]
The only positive roots are
\[
r_1 = \left( \frac{4 - 2 \sqrt{4 - 5e}}{e} \right)^{1/6} \;, \qquad r_2 = \left( \frac{4 + 2 \sqrt{4 - 5e}}{e} \right)^{1/6}
\]

Actually we know that only one root is correct. The last condition tells us that
\[
f''_ c(r) + \frac 3 r f'_ c(r) = \frac{96}{e \, r^{14}} (-5 + r^6) > 0
\]
so $r > 5^{1/6}$. So $r_1$ is unqualified for this condition:
\[
\begin{align}
r_1^6 &= \frac{4 - 2 \sqrt{4 - 5e}}{e} = \frac{4 - 2 \sqrt{4 - 5e}}{e} \cdot \frac{4 + 2 \sqrt{4 - 5e}}{4 + 2 \sqrt{4 - 5e}} \\
&= \frac{20e}{e (4 + 2 \sqrt{4 - 5e})} < \frac{20}{4+0} = 5
\end{align}
\]
By the same way, we know that $r=r_2$ is the only root qualified for the last condition provided.

We then substitute $r=r_2$ into the first condition to give $b$. The only positive root expression of $b$ is something awful:
\[
b_\mathrm{crit} = \frac{2^{2/3} 3^{1/2}}{5} \left( \frac{2 + \sqrt{4-5e}}{e} \right)^{1/6} \left( \frac{-2+ \sqrt{4-5e}-5e}{e} \right)^{1/2}
\]

From the solutions mentioned above, only when
\[
0 < e < \frac 4 5
\]
can the central particle "trap" be formed.
