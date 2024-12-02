# Project - FEM
## B.1 Heat Equilibrium

In this project, we wish to simulate heat conduction in a non-insulated water hose filled with water that is assumed to be at rest. The surrounding of the hose has a temperature of \( u_0 \, \text{K} \). We assume that the water in the hose is hit by microwaves that cause the water to heat up. The microwaves act as heat sources because their energy turns into heat when the waves are absorbed by the water. This heat source is described by a function \( f(x) \). 

We consider a cross-section of the water hose and consider the stationary heat equation. Let the hose’s cross-section be described by a domain \( D \) and let \( u : D \to \mathbb{R} \) be the heat distribution. Then the temperature distribution is described by the boundary value problem:

\[
\begin{cases} 
-\text{div}(a \nabla u(x)) = f(x), & x \in D, \\
a \partial_\nu u(x) = c (u_0 - u(x)), & x \in \partial D.
\end{cases} \tag{B.1}
\]

Here \( a \) is the thermal conductivity coefficient of water, and \( c \) is the heat conductivity coefficient of the hose’s walls. Typical values are:

- \( a_{\text{water}} = 0.6 \, \text{W/K} \),
- \( c = 10 \, \text{W/(dm K)} \).

We assume that the cross-section of the hose is a circle of radius 1 dm. Note that we use dm as the unit of length, so \( r = 1 \) and not \( r = 0.1 \)!

---

### Tasks

#### 1. Variational Formulation and FEM Implementation
- **Task:** Find the variational formulation \( V \), \( a \), and \( L \) of the boundary value problem and write a MATLAB program that solves the BVP using FEM. The finite element integrals must be computed by hand. Take for now \( a, c, f, u_0 \) as constants, but allow \( D \) to be an arbitrary domain (although we will always use the above disk \( D \)).
- **Report:** Pseudocode for functions `IntMatrix`, `BdyMatrix`, `IntVector` and `BdyVector`, and computations of finite element integrals.

---

#### 2. Code Validation Against Analytic Solution
- **Task:** We next check the code against an analytic solution. Consider the function:

\[
u(x, y) = 325 - 20(x^2 + y^2),
\]

where \( a = 0.6 \), \( c = 10 \), \( u_0 = 300 \). Find the corresponding values for \( c \) and \( f \). Use your code to calculate \( u \) numerically for these \( a, c, f, u_0 \). Check your code if the error does not seem to converge to 0.  
*(Hint: Carefully specify your circle with `pdecirc`.)*

- **Report:** Maximum absolute errors of the solution for different meshes. The rate of convergence \( O(h^N) \): What is the integer \( N \)?

---

#### 3. Non-linear PDE Problem
- **Task:** Consider the following non-linear PDE problem. Set the outside temperature to \( u_0 = 250 \) and \( c = 10 \). Ice has a different thermal conductivity coefficient than water, namely \( a_{\text{ice}} = 2.2 \, \text{W/K} \). So \( a = a(u) \) depends on the temperature \( u \) now:

\[
a(u) = 
\begin{cases} 
2.2, & u < 273, \\
0.6, & u > 273.
\end{cases} \tag{B.2}
\]

This gives a non-linear partial differential equation that we can solve with a fixed-point iteration.  
*(See Section 2.5 for details on fixed-point iteration.)*

The procedure is:
1. Assume first \( a = 0.6 \) everywhere.
2. Solve the boundary value problem and adjust \( a \) according to the values of \( u \). 
3. Repeat this process as long as \( a \) has to be adjusted.

What constant value of the heat source function \( f \) is needed so that in the equilibrium state, we have 50% ice and 50% water of the area cross-section?

- **Report:** The value of \( f \) and a plot of the non-linear solution \( u \).

---

### Extra Exercises (1 Bonus Point Each)

#### 4. Non-Constant \( f \) Using Quadrature (D.2)
(a) We check the new code using the function:

\[
u(x, y) = 325 - 20(x^2 + y^2)^2,
\]

where \( a = 0.6 \), \( u_0 = 300 \). Find the corresponding values for \( c \) and \( f \). Use your code to calculate \( u \) numerically and compare the result to the above exact solution. Check your code if the error does not seem to converge to 0.

(b) Let \( a = 0.6 \), \( c = 10 \), \( u_0 = 300 \), and use the source:

\[
f(x, y) = f_0 e^{-\mu \sqrt{1 - y^2 - x}},
\]

arising from the Beer-Lambert law with attenuation \( \mu = \ln 2 \) and source incident from the right. Find the value of \( f_0 \) that gives a maximal temperature of 320 in the hose.

- **Report:** 
  - Pseudocode for the assembly with `IntVector`.
  - (a) Error estimates for different meshes.
  - (b) The value of \( f_0 \) and a plot of the corresponding heat distribution.

---

#### 5. Log-Plot Analysis of Condition Number
- **Task:** Let \( a = 0.6 \), \( u_0 = 300 \), and \( f(x, y) = 100 \sin(6y) \). Vary the value of \( c \):
  - (a) Log-plot the condition number of \( A \) for your discretized system \( Ax = b \), as a function of \( -\infty < c < +\infty \). Find more than one value of \( c \) for which \( A \) is not invertible/ill-conditioned.
  - (b) When is your bilinear form \( a(u, \phi) \) coercive?
  - (c) Plot the solution \( u \) when:
    - \( c \approx 0 \),
    - \( c \approx +\infty \),  
    and explain the results.

- **Report:** The plots and discussion.



## Project - BIE
TBC
