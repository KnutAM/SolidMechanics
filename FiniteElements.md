+++
 title = "Finite elements"
+++

**Intended learning outcomes**
* Derive the weak form of the equilibrium from the strong form of the equilibrium equations
* Derive the finite element form from the weak form
* Calculate strains, stresses and element nodal force contributions for simple elements.

# The Finite Element Method (FEM)
We use the finite element method to solve partial differential equations describing a boundary value problem. 
Here, we discuss the solution of the [mechanical boundary value problem](/BoundaryValueProblem). 
The goal is only to briefly introduce the method, and not to discuss the many interesting aspects related to its implementation 
and use in various settings. This is left for specialized courses. 

## Strong form
The strong form of our boundary value problem is described [here](/BoundaryValueProblem), and summarized by the following figure and equations.

![](/assets/Potato.svg)

\begin{align}
\div{\sig(\eps(\tv{u}(\tv{x})))} + \tv{b} &= \tv{0},\phantom{(\tv{x})} \quad \tv{x} \in \Omega \\
\tv{u}(x) &= \tv{u}\subscr{D}(\tv{x}), \quad \tv{x} \in \Gamma\subscr{D} \\
\tv{n} \cdot \sig(\eps(\tv{u}(\tv{x}))) &= \tv{t}\subscr{N}(\tv{x}), \quad \tv{x} \in \Gamma\subscr{N}\label{eq:strongform}
\end{align}

## Weak form
The weak form is an integrated form of the strong form, which makes the problem possible to approximate using the [Finite Element Method](#finite_element_form). To obtain this form, we apply two "tricks":

1. Scalar multiply the strong form with an arbitrary test displacement field $\tv{\delta u}(\tv{x})$
2. Integrate over the domain $\Omega$

Because the test displacement field $\tv{\delta u}(\tv{x})$ is arbitrary, the following weak form of the boundary value problem is in practice equivalent to Equation \eqref{eq:strongform}. (There are some restrictions wrt. to the smoothness of the functions $\tv{u}(\tv{x})$ and $\tv{\delta u}(\tv{x})$ from mathematics.) For brevity, from hereon we don't write out the dependence on $\tv{x}$ and $\eps$ for $\tv{u}$, $\tv{\delta u}$, $\sig$, $\tv{u}\subscr{D}$, and $\tv{t}\subscr{N}$. 
\begin{align}
\int_\Omega \tv{\delta u} \cdot \left[\div{\sig}\right]\, \dif \Omega &= - \int_\Omega \tv{\delta u} \cdot \tv{b}\, \dif \Omega \\
\tv{u} &= \tv{u}\subscr{D} \; \mathrm{on}\, \Gamma\subscr{D} \\
\tv{n} \cdot \sig(\eps(\tv{u})) &= \tv{t}\subscr{N} \; \mathrm{on}\, \Gamma\subscr{N}
\end{align}

If we then apply the [Green-Gauss theorem](https://knutam.github.io/tensors/Theory/TensorIntegration/) to the first equation, we obtain
\begin{align}
\int_\Gamma \tv{\delta u}\cdot \sig \cdot \tv{n}\, \dif \Gamma - \int_\Omega \left[\grad{\tv{\delta u}}\right] : \sig\, \dif \Omega &= - \int_\Omega \tv{\delta u} \cdot \tv{b}\, \dif \Omega \\
\end{align}
and we can insert the known traction on the Neumann boundary $\Gamma\subscr{N}$ as
\begin{align}
\int_{\Gamma\subscr{N}} \tv{\delta u}\cdot \tv{t}\subscr{N}\, \dif \Gamma\subscr{N} + \int_{\Gamma\subscr{D}} \tv{\delta u}\cdot \sig \cdot \tv{n}\, \dif \Gamma\subscr{D} - \int_\Omega \left[\grad{\tv{\delta u}}\right] : \sig\, \dif \Omega &= - \int_\Omega \tv{\delta u} \cdot \tv{b}\, \dif \Omega \\
\tv{u} &= \tv{u}\subscr{D} \; \mathrm{on}\, \Gamma\subscr{D} \label{eq:weakform}
\end{align}
while keeping the known displacements on $\Gamma\subscr{D}$. 

## Finite element form
The weak form in Equation \eqref{eq:weakform} is also very difficult to solve, perhaps arguably more difficult than the strong form in Equation \eqref{eq:strongform} since we need to deal with the arbitrary function $\tv{\delta u}$. The main issue for both of these forms, is that the solution is a continuous function $\tv{u}$. It would be much easier to solve if we could somehow parameterize the function, e.g. if we could write that (in 1D) $u(x) = a_1 x + a_2 x^2 + a_3 x^3$ and then just determine the coefficients $a_i$. This is what we will do, namely that we will approximate both $\tv{u}$ and $\tv{\delta u}$ using so-called base functions:
\begin{align}
\tv{u}(\tv{x}) &\approx \tv{N}_i(\tv{x}) a_i \\
\tv{\delta u}(\tv{x}) &\approx \tv{N}^\delta_i(\tv{x}) \delta a_i
\end{align}
where $a_i$ and $\delta a_i$ are the coefficients that we need to determine and $\tv{N}_i$ and $\tv{N}^\delta_i$ are shape functions that we postulate can be used to describe the solution and test space. So we restrict our solution to solutions that can be described by this *discretization*. Here, we will use the standard Galerkin method of having $\tv{N}_i = \tv{N}^\delta_i$, which is described in detail in specific finite element courses. Inserting these approximations into our weak form, we obtained the discretized FE-problem 
\begin{align}
\delta a_i \left[\int_{\Gamma\subscr{N}} \tv{N}_i \cdot \tv{t}\subscr{N}\, \dif \Gamma\subscr{N} + \int_{\Gamma\subscr{D}}\tv{N}_i \cdot \sig \cdot \tv{n}\, \dif \Gamma\subscr{D} - \int_\Omega \left[\grad{ \tv{N}_i}\right] : \sig\, \dif \Omega\right] &= - \delta a_i \int_\Omega \tv{N}_i \cdot \tv{b}\, \dif \Omega \\
a_i \tv{N}_i = \tv{u}\subscr{D} \; \mathrm{on}\, \Gamma\subscr{D} \label{eq:feform}
\end{align}
Since $\tv{\delta u}$ should be arbitrary, so should $\delta a_i$. So by letting $\delta a_i=1$ and $\delta a_j=0$ for all $j\neq i$ for all $i$ (e.g. $\delta a_1=1$, then $\delta a_2=\delta a_3=\cdots=\delta a_n=0$ and $\delta a_2=1$, then $\delta a_1=\delta a_3=\cdots=\delta a_n=0$ and so on), we have $n$ equations ($i=1,2,\cdots,n$) that we can write as
\begin{align}
f\supscr{int}_i &= f\supscr{ext}_i \\
f\supscr{int}_i &= \int_\Omega \left[\grad{ \tv{N}_i}\right] : \sig\, \dif \Omega \\
f\supscr{ext}_i &= \int_{\Gamma\subscr{N}} \tv{N}_i \cdot \tv{t}\subscr{N}\, \dif \Gamma\subscr{N} + \int_{\Gamma\subscr{D}}\tv{N}_i \cdot \sig \cdot \tv{n}\, \dif \Gamma\subscr{D} + \int_\Omega \tv{N}_i \cdot \tv{b}\, \dif \Omega \label{eq:feform}
\end{align}

where $f\supscr{int}$ and $f\supscr{ext}$ are the so-called internal (due to stresses) and external (due to applied loads) forces, respectively. 

### Base functions
The base functions $\tv{N}$ have many names, e.g. *basis*, *shape*, *test*, *trial*, etc. functions. These are generally used interchangably, potentially leading to some confusion.
In general, *test* and *trial* often denote the functions approximating the *test* function, $\tv{\delta u}$ (yes, the same name), i.e. $\tv{N}^\delta$. *Shape*
functions often describe the functions approximating the solution $\tv{u}$, i.e. $\tv{N}$. Finally, *base* or *basis* functions are often used to describe both $\tv{N}$ and $\tv{N}^\delta$.
For the $n$ equations $f\supscr{int}_i = f\supscr{ext}_i$ above to form a solvable equation system, we require (1) that the base functions are independent. That is, for each $j\in[1,n]$, there doesn't exist 
coefficients $a_i$ where $a_j=0$ such that $\tv{N}_j = a_i \tv{N}_i$. Secondly, we require that the solution is constrained by Dirichlet boundary conditions such that rigid body motions are prevented. 

In practice, we often choose so-called nodal base functions. Consider the nodal positions $\tv{x}_i$. We then have a set of scalar base functions $N_i$, which are defined such that
\begin{align}
N_i(\tv{x}_j) = \left\lbrace \begin{matrix} 1 & i=j \\ 0 & i\neq j \end{matrix}\right.
\end{align}

In one dimension for linear base functions, this choice looks like

![](/assets/BaseFunctions1D.svg#halfwidth)

for 3 elements. If we consider the definition of the base function inside one element, between nodes $x_i$ and $x_{i+1}$, we have
\begin{align}
N_i(x) = 1 - \frac{x-x_i}{x_{i+1}-x_{i}}, \quad N_{i+1}(x) = \frac{x-x_i}{x_{i+1}-x_{i}}
\end{align}

As we increase to two dimensions, we have a more complicated situation. The left figure below shows only the triangular element, and then each of the base function in the other three figures.

![](/assets/BaseFunctions2D.svg#fullwidth)

The linear base functions inside an element with nodal coordinates $\tv{x}_1$, $\tv{x}_2$, and $\tv{x}_3$, ordered counter-clockwise around the element as in the figure above, are
\begin{align}
N_1(x,y) &= \frac{x_2 y_3 - x_3 y_2 + [y_2-y_3]x + [x_3-x_2]y}{2A} \\
N_2(x,y) &= \frac{x_3 y_1 - x_1 y_3 + [y_3-y_1]x + [x_1-x_3]y}{2A} \\
N_3(x,y) &= \frac{x_1 y_2 - x_2 y_1 + [y_1-y_2]x + [x_2-x_1]y}{2A}
\end{align}
where $A=\left[[x_2 y_3 - x_3 y_2] - [x_1 y_3 - x_3 y_1] + [x_1 y_2 - x_2 y_1]\right]/2$ is the element's area. It is left as an exercise for the reader to verify that $N_i(x_j,y_j)=\delta_{ij}$. For a detailed derivation, please see *Niels Ottosen and Hans Petersson, 1992: "Introduction to the Finite Element Method"*, Section 7.3.1. 

However, the base functions discussed so far are scalar base functions. But since we want to approximate vector functions, we need vector base functions. These are easily constructed for each node as
\begin{align}
\tv{N}_{[3i-2]}(\tv{x}) &= [N_i(\tv{x}), 0, 0]\\
\tv{N}_{[3i-1]}(\tv{x}) &= [0, N_i(\tv{x}), 0] \\
\tv{N}_{[3i]}(\tv{x}) &= [0, 0, N_i(\tv{x})]
\end{align}
for three dimensions. More formally and short, we can write $\tv{N}_i(\tv{x}) = N_j(\tv{x})\onebase{k}$ where $j$ is the node number, $m$ is the dimension, $1\leq k\leq m$, and $i=m(j-1)+k$. If we have $n\subscr{nod}$ nodes, we then have $n=m\,n\subscr{nod}$ base functions $\tv{N}_i$. And as mentioned above, there are equally many coefficients $a_i$ and equations $\tv{f}\supscr{int}_i=\tv{f}\supscr{ext}_i$ as base functions.


### Evaluating the integrals
The force vectors $f\supscr{int}_i$ and $f\supscr{ext}_i$ in Equation \eqref{eq:feform} are calculated by solving the integral over the entire domain. But the construction of base functions above ensures that one base function is only non-zero in a few elements, as illustrated below. 

![](/assets/BaseFunctionsMesh.svg#fullwidth)

This makes it possible to calculate contributions elementwise. If we evaluate $\grad{\tv{N}_i}$, for the base functions above, this becomes a constant tensor as $N_j(\tv{x})$ is linear inside the element. Consequently, the strain $\tv{\epsilon} = \left[\grad{\tv{u}}\right]\sym = a_i \left[\grad{\tv{N}_i}\right]\sym$ is also constant. And since the stress is a function of the strain, then $\sig(\eps)$ is constant inside the element, implying that the entire integrand $\left[\grad{ \tv{N}_i}\right] : \sig$ is constant. So for triangular elements with linear base functions, the integrals are trivial to solve. For cases when the base functions are not linear, e.g. quadrilateral elements (4 nodes) or quadratic (6-noded) triangular elements, we need to use numerical integration techniques, see e.g. *Niels Ottosen and Hans Petersson, 1992: "Introduction to the Finite Element Method"*, Chapter 20. 