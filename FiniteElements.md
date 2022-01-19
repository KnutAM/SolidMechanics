+++
 title = "Finite elements"
+++

**Intended learning outcomes**
* Derive the weak form of the equilibrium from both the principle of minimum potential energy and the strong form of the equilibrium equations
* Derive the finite element form from the weak form

# Finite elements
Intro...

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

1. Scalar multiply the strong form with an arbitrary test displacement field $\tv{v}(\tv{x})$
2. Integrate over the domain $\Omega$

Because the test displacement field $\tv{v}(\tv{x})$ is arbitrary, the following weak form of the boundary value problem is in practice equivalent to Equation \eqref{eq:strongform}. (There are some restrictions wrt. to the smoothness of the functions $\tv{u}(\tv{x})$ and $\tv{v}(\tv{x})$ from mathematics.) For brevity, from hereon we don't write out the dependence on $\tv{x}$ and $\eps$ for $\tv{u}$, $\tv{v}$, $\sig$, $\tv{u}\subscr{D}$, and $\tv{t}\subscr{N}$. 
\begin{align}
\int_\Omega \tv{v} \cdot \left[\div{\sig}\right]\, \dif \Omega &= - \int_\Omega \tv{v} \cdot \tv{b}\, \dif \Omega \\
\tv{u} &= \tv{u}\subscr{D} \; \mathrm{on}\, \Gamma\subscr{D} \\
\tv{n} \cdot \sig(\eps(\tv{u})) &= \tv{t}\subscr{N} \; \mathrm{on}\, \Gamma\subscr{N}
\end{align}

If we then apply the [Green-Gauss theorem](https://knutam.github.io/tensors/Theory/TensorIntegration/) to the first equation, we obtain
\begin{align}
\int_\Gamma \tv{v}\cdot \sig\trans \cdot \tv{n}\, \dif \Gamma - \int_\Omega \left[\grad{\tv{v}}\right] : \sig\trans\, \dif \Omega &= - \int_\Omega \tv{v} \cdot \tv{b}\, \dif \Omega \\
\end{align}
and we can insert the known traction on the Neumann boundary $\Gamma\subscr{N}$ as
\begin{align}
\int_{\Gamma\subscr{N}} \tv{v}\cdot \tv{t}\subscr{N}\, \dif \Gamma\subscr{N} + \int_{\Gamma\subscr{D}} \tv{v}\cdot \sig\trans \cdot \tv{n}\, \dif \Gamma\subscr{D} - \int_\Omega \left[\grad{\tv{v}}\right] : \sig\trans\, \dif \Omega &= - \int_\Omega \tv{v} \cdot \tv{b}\, \dif \Omega \\
\tv{u} &= \tv{u}\subscr{D} \; \mathrm{on}\, \Gamma\subscr{D} \label{eq:weakform}
\end{align}
while keeping the known displacements on $\Gamma\subscr{D}$. 

## Finite element form
The weak form in Equation \eqref{eq:weakform} is also very difficult to solve, perhaps arguably more difficult than the strong form in Equation \eqref{eq:strongform} since we need to deal with the arbitrary function $\tv{v}$. The main issue for both of these forms, is that the solution is a continuous function $\tv{u}$. It would be much easier to solve if we could somehow parameterize the function, e.g. if we could write that (in 1D) $u(x) = a_1 x + a_2 x^2 + a_3 x^3$ and then just determine the coefficients $a_i$. This is what we will do, namely that we will approximate both $\tv{u}$ and $\tv{v}$ using so-called shape functions:
\begin{align}
\tv{u}(\tv{x}) &\approx \tv{N}\supscr{u}_i(\tv{x}) a_i
\tv{v}(\tv{x}) &\approx \tv{N}\supscr{v}_i(\tv{x}) \delta a_i
\end{align}
where $a_i$ and $\delta a_i$ are the coefficients that we need to determine and $\tv{N}\supscr{u}_i$ and $\tv{N}\supscr{v}_i$ are shape functions that we postulate can be used to describe the solution and test space. So we restrict our solution to solutions that can be described by this *discretization*. Here, we will use the standard Galerkin method of having $\tv{N}_i = \tv{N}\supscr{u}_i = \tv{N}\supscr{v}_i$, which is described in detail in specific finite element courses. Inserting these approximations into our weak form, we obtained the discretized FE-problem 
\begin{align}
\delta a_i \left[\int_{\Gamma\subscr{N}} \tv{N}_i \cdot \tv{t}\subscr{N}\, \dif \Gamma\subscr{N} + \int_{\Gamma\subscr{D}}\tv{N}_i \cdot \sig\trans \cdot \tv{n}\, \dif \Gamma\subscr{D} - \int_\Omega \left[\grad{ \tv{N}_i}\right] : \sig\trans\, \dif \Omega\right] &= - \delta a_i \int_\Omega \tv{N}_i \cdot \tv{b}\, \dif \Omega \\
a_i \tv{N}_i = \tv{u}\subscr{D} \; \mathrm{on}\, \Gamma\subscr{D} \label{eq:feform}
\end{align}
Since $\tv{v}$ should be arbitrary, so should $\delta a_i$. So by letting $\delta a_i=1$ and $\delta a_j=0$ for all $j\neq i$ for all $i$ (e.g. $\delta a_1=1$, then $\delta a_2=\delta a_3=\cdots=\delta a_n=0$ and $\delta a_2=1$, then $\delta a_1=\delta a_3=\cdots=\delta a_n=0$ and so on), we have $n$ equations ($i=1,2,\cdots,n$) that we can write as
\begin{align}
f\supscr{int}_i &= f\supscr{ext}_i \\
f\supscr{int}_i &= \int_\Omega \left[\grad{ \tv{N}_i}\right] : \sig\trans\, \dif \Omega \\
f\supscr{ext}_i &= \int_{\Gamma\subscr{N}} \tv{N}_i \cdot \tv{t}\subscr{N}\, \dif \Gamma\subscr{N} + \int_{\Gamma\subscr{D}}\tv{N}_i \cdot \sig\trans \cdot \tv{n}\, \dif \Gamma\subscr{D} + \int_\Omega \tv{N}_i \cdot \tv{b}\, \dif \Omega 
\end{align}

where $f\supscr{int}$ and $f\supscr{ext}$ are the so-called internal (due to stresses) and external (due to applied loads) forces, respectively. 

## Implementation
To be completed...