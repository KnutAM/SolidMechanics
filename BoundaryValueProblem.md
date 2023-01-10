+++
 title = "Boundary value problem"
+++

**Intended learning outcomes**
* Derive the strong form of the equilibrium equations for a solid
* Apply appropriate boundary conditions to the boundary value problem


# The boundary value problem (BVP)
Our goal is to investigate a structure subjected to external loads. We might, for example, be interested in the resulting deformations or the internal stresses. One way to think about this, is that applying known displacements to the surface causes stresses, $\sig$, inside the body. For equilibrium to hold, we require that these stresses fulfill
\begin{align}
\div{\sig} + \tv{b} &= \tv{0}
\end{align}
inside the body subjected to the body load $\tv{b}$ [Force/Volume], and that the traction $\tv{t}$ is
\begin{align}
\tv{t} &= \tv{n}\cdot\sig
\end{align}
on the boundary, where $\tv{n}$ is the normal vector. 

![](/assets/Potato.svg)

How high the stresses become for given displacements, depends on the material behavior. That is, how $\sig$ depends on the strain $\eps$. For linear elasticity, we have that
\begin{align}
\sig(\eps) = \tf{E}:\eps, \quad \text{(Linear elasticity)}
\end{align}

We thus need to determine the strains, $\eps$, inside the body to calculate the stresses. And the strains depend on the displacements $\tv{u}$ as
\begin{align}
\eps(\tv{u}) = \left[\grad{\tv{u}}\right]\sym
\end{align}

We now have the differential equations required to solve the boundary value problem, namely 
\begin{align}
\div{\sig(\eps(\tv{u}(\tv{x})))} + \tv{b} &= \tv{0},\phantom{(\tv{x})} \quad \tv{x} \in \Omega \\
\tv{u}(x) &= \tv{u}\subscr{D}(\tv{x}), \quad \tv{x} \in \Gamma\subscr{D} \\
\tv{n} \cdot \sig(\eps(\tv{u}(\tv{x}))) &= \tv{t}\subscr{N}(\tv{x}), \quad \tv{x} \in \Gamma\subscr{N}\label{eq:bvp}
\end{align}
where the first equation defines the differential equation inside the body (light orange). The second equation give the so-called Dirichlet boundary conditions (black, $\Gamma\subscr{d}$), here we already know the value of the displacements. And finally, the part of the boundary with unknown displacements is called the Neuman (also natural) boundary (blue, $\Gamma\subscr{N}$). Here, we know the traction. Often, a majority of the boundary has zero load but we know it is zero!

The problem is - it is not possible to determine the function $\tv{u}(\tv{x})$ for most cases. While it is possible in a few special cases, what to do for the remaining problems? Often, we use [the finite element method](/FiniteElements). But now let's first start with a case where we can find the solution analytically.

## Analytical solution to the BVP (uniaxial stress)
Consider the case of cylinder subjected to uniaxial tension, with an isotropic material:

![](/assets/UniaxialTensionTest.svg)

On the blue edges we have Neumann (natural) boundary conditions. We don't know the displacements, but we know the traction. On the right side ($x_1=L$), this traction is $[F/A]\onebase{1}$, where $A$ is the cross-sectional area and on the sides of the cylinder the traction is zero. On the left side, we actually have mixed boundary conditions. We know that the traction components in the $x_2$-$x_3$ plane are zero (Neumann), and that the displacement is zero in the $\onebase{1}$ direction (Dirichlet). 

\collaps{To solve the differential equation, we make the following ansatz}{**How did we arrive at this ansatz?**

We assume that the deformation will be homogeneous, i.e. there will be no shear and hence each "disk" with normal $\onebase{1}$ will deform equally. Firstly, we know that the isotropic material will not shear if loaded only in a normal direction. Secondly, we see that both sides of the cylinder behave identically with only loads in the normal direction and no shear traction forces).}
\begin{align}
\tv{u}(\tv{x}) = a_1 x_1 \onebase{1} + a_2 x_2 \onebase{2} + a_3 x_3 \onebase{3}
\end{align}

From this ansatz, we obtain the constant strains as
\begin{align}
\eps = a_1 \twobase{1}{1} + a_2 \twobase{2}{2} + a_3 \twobase{3}{3} \label{eq:strainansatz}
\end{align}

For linear isotropic elasticity, we have $\sig = 2G \eps\dev + 3K\eps\sph$, with
\begin{align}
\eps\sph &= \frac{a_1+a_2+a_3}{3} \ts{I} \\
\eps\dev &= \frac{2a_1 - a_2 - a_3}{3}\twobase{1}{1} + \frac{2a_2 - a_1 - a_3}{3}\twobase{2}{2} + \frac{2a_3 - a_1 - a_2}{3}\twobase{3}{3}
\end{align}

Checking the equilibrium by taking the divergence of $\sig$, we see that this is zero as $\sig$ is constant. Hence, the proposed ansatz fulfills the equilbrium equation. However, we must also fulfill the boundary conditions. At $x_1=L$, we have that $\tv{t}=[F/A]\onebase{1}$, that is
\begin{align}
\tv{t} &= \onebase{1} \cdot \sig = \onebase{1} \cdot [2G \eps\dev + 3K\eps\sph] \\
&= 2G \frac{2a_1 - a_2 - a_3}{3} \onebase{1} + K [a_1+a_2+a_3] \onebase{1} = \frac{F}{A}\onebase{1} \\
\frac{F}{A} &= 2G \frac{2a_1 - a_2 - a_3}{3} + K [a_1+a_2+a_3]
\end{align}

Considering the boundary conditions on the side of the cylinder, e.g. at $x_2=r,\,x_3=0$ and $x_2=0,\,x_3=r$, we have that the traction is zero, e.g. 
\begin{align}
\tv{t} &= \onebase{2} \cdot \sig = 2G \frac{2a_2 - a_1 - a_3}{3} \onebase{2} + K [a_1+a_2+a_3] \onebase{2} = \tv{0} \\
0 &= 2G \frac{2a_2 - a_1 - a_3}{3} + K [a_1+a_2+a_3] \\
0 &= 2G \frac{2a_3 - a_1 - a_2}{3} + K [a_1+a_2+a_3], \quad \mathrm{obtained\,the\,same\,way\,as\,above}
\end{align}

And hence, we see from the two last equivalent equations where we just swapped $a_2$ and $a_3$ that $a_2=a_3$. Resulting in 
\begin{align}
\frac{F}{A} &= 4G \frac{a_1 - a_2}{3} + K [a_1+2a_2] \\
0 &= 2G \frac{a_2 - a_1}{3} + K [a_1+2a_2] \\
a_2 &= -\frac{3K - 2G}{2G + 6K}a_1
\end{align}
To simplify, further calculation, we denote $\nu = [3K-2G]/[2G+6K]$ (This is the actual Poissons' ratio, as we can see from how it governs the relation between $a_1$ and $a_2$ along with Equation \eqref{eq:strainansatz}). We then have
\begin{align}
a_2 &= -\nu a_1 \\
\frac{F}{A} &= 4G \frac{a_1 + \nu a_1}{3} + K [a_1-2\nu a_1] \\
a_1 &= \frac{F}{A[4G[1+\nu] + K[1-2\nu]]}
\end{align}
This expression is now simplified by denoting $E=4G[1+\nu] + K[1-2\nu]$ (This is the actual Young's modulus, as we will se later). And we have that
\begin{align}
a_1 &= \frac{F}{AE} \\
a_2 &= -\nu \frac{F}{AE} = a_3
\end{align}
Yielding the following displacements and strains. 
\begin{align}
\tv{u}(\tv{x}) &= \frac{F}{AE} \left[ x_1 \onebase{1} - \nu x_2 \onebase{2} - \nu x_3 \onebase{3}\right] \\
\eps &= \frac{F}{AE} \left[ \twobase{1}{1} - \nu [\twobase{2}{2} + \twobase{3}{3}]\right] \label{eq:kinematicresult}
\end{align}
Finally, if we insert the constitutive relationship $\sig=2G \eps\dev + K\eps\sph$, and using our definitions of $E$ and $\nu$ from above, we obtain
\begin{align}
\sig &= \frac{F}{A} \twobase{1}{1} \label{eq:stressresult}
\end{align}
I.e. we see that the stress is indeed uniaxial. Furthermore, seeing combining Equation \eqref{eq:kinematicresult} with Equation \eqref{eq:stressresult}, we see that $\sigma_{11}=E \epsilon_{11}$, showing that the $E$ we defined earlier is the actual Young's modulus. 

In summary, for even this very simple case, determining the solution to the BVP in Equation \eqref{eq:bvp} was quite a bit of effort. And as previously noted, in many cases no unique solution exists even if we include simplifications such as plane stress or strain.
