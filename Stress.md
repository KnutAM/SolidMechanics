+++
 title = "Stress"
+++

**Intended learning outcomes**
* Define the traction vector
* Explain Cauchy's theorem
* Define the Cauchy stress
* Derive the continuum force equilibrium equation
* Show that the Cauchy stress is symmetric

# Stress
It is difficult to give a short answer to the question "what is stress?". Here, we will explain it via the more understandable traction vector. Then, via Cauchy's theorem, we can introduce the notion Cauchy stress tensor $\sig$. With the definition of $\sig$, we can derive the equilibrium equation for a continuum body. Finally, we will show that equilibrium also requires that $\sig$ is a symmetric tensor. 

## The traction vector
![](/assets/traction.svg)

Consider that we have a loaded structure (left), and then we make a plane cut through this structure (green line). In order to replace the internal forces acting on the cut plane, we need to introduce forces. These can be both perpendicular to the cut plane and along the cut plane. In reality, the forces are discrete between atoms in the material. However, in continuum mechanics, we introduce this as a continuous field of force per area, which we define as traction. More precisely, we consider an infitesimal area $\dif A$ of this cut surface with normal vector $\tv{n}$ (right). The sum of force vectors acting on this surface is $\tv{\dif F}$ and the traction is then defined as
\begin{align}
\tv{t} := \frac{\tv{\dif F}}{\dif A}
\end{align}



## Cauchy's theorem
Above, we defined the traction vector $\tv{t}$. For a loaded body, the traction vector, $\tv{t}$, depends on both the position in the body and the normal vector of the cut plane. It would be convenient to have a measure of the internal forces (per area) in the body that doesn't depend on the exact cut plane. To solve this, we consider the tetrahedron below. 

![](/assets/StressTetrahedron.svg)

### Relationship between areas
To be able to use this, we need to determine the relationship between the four areas. To this end, we will use the divergence theorem. If we consider a vector function $\tv{a}(\tv{x})$ defined inside and in the viscinity of a closed body $\Omega$ with the boundary $\Gamma$, the Gauss (divergence) theorem states
\begin{align}
\int_\Omega \div{\tv{a}} \dif \Omega &= \int_\Gamma \tv{a}\cdot\tv{n} \dif \Gamma \\
\int_\Omega a_j \nabla_j \dif \Omega &= \int_\Gamma a_j n_j \dif \Gamma
\end{align}
where $\tv{n}$ is a unit normal vector pointing out of the body. By using this equation three times with $b_{1j}=a_j$, $b_{2j}=a_j$, and lastly $b_{3j}=a_j$ we get three equations as
\begin{align*}
\int_\Omega b_{ij} \nabla_j \dif \Omega &= \int_\Gamma b_{ij} n_j \dif \Gamma
\end{align*}
If we then set $b_{ij}=\delta_{ij}$, we obtain
\begin{align}
\int_\Gamma n_i \dif \Gamma = 0
\end{align}
as $\delta_{ij} \nabla_j$ is zero as $\delta_{ij}$ is constant. 
For our tetrahedron, this implies that
\begin{align}
A \tv{n} = A_1 \onebase{1} + A_2 \onebase{2} + A_3 \onebase{3} = A_i \basei
\end{align}
where $\tv{n}$ is the normal vector of the plane with area $A$. The three other faces have normal vectors $\tv{n}_1=-\onebase{1}$, $\tv{n}_2=-\onebase{2}$, and $\tv{n}_3=-\onebase{3}$. Scalar multiplying this expression by $\onebase{j}$ we obtain 
\begin{align}
A \tv{n}\cdot\onebase{j} = A_i \delta_{ij} = A_j \label{eq:arearelation}
\end{align}

### Equilibrium on tetrahedron
We consider that the tetrahedron can include a body load which is load per volume, $\tv{b}$, where the volume is $V$. Then, the equilibrium equation for the tetrahedron becomes
\begin{align*}
0 = \tv{t}_i A_i + \tv{t} A + \tv{b} V
\end{align*}

Evaluating this equation for each coordinate direction, $\onebase{j}$, we get
\begin{align*}
0 &= A_i \tv{t}_i \cdot \onebase{j} + A \tv{t} \cdot \onebase{j} + V \tv{b} \cdot \onebase{j}
\end{align*}
Inserting Equation \eqref{eq:arearelation}, we obtain
\begin{align*}
0 &= A [\tv{n}\cdot\onebase{i}] [\tv{t}_i \cdot \onebase{j}] + A \tv{t} \cdot \onebase{j} + V \tv{b} \cdot \onebase{j} \\
0 &= [\tv{n}\cdot\onebase{i}] [\tv{t}_i \cdot \onebase{j}] + \tv{t} \cdot \onebase{j} + \frac{V}{A} \tv{b} \cdot \onebase{j}
\end{align*}
If we now let the tetrahedron shrink, the volume to area ratio goes to zero (because volume is proportional to the side lengths cubed and area to lengths squared). Hence, we get
\begin{align*}
[\tv{n}\cdot\onebase{i}] [\tv{t}_i \cdot \onebase{j}] + \tv{t} \cdot \onebase{j} = 0
\end{align*}
Let us denote the quantity $-\tv{t}_i \cdot \onebase{j}=\sigma_{ij}$. Noting that $n_i=\tv{n}\cdot\onebase{i}$ is the component of $\tv{n}$ in the $\onebase{i}$ direction, and $t_j=\tv{t} \cdot \onebase{j}$ is the component of $\tv{t}$ in the $\onebase{j}$ direction, we have 
\begin{align*}
n_i \sigma_{ij} = t_j
\end{align*}
Or, in tensor form
\begin{align}
\tv{n}\cdot \sig = \tv{t} \label{eq:CauchysTheorem}
\end{align}
Equation \eqref{eq:CauchysTheorem} is **Cauchy's Theorem**. We see that the 2nd order tensor $\sig$ can describe the traction on plane at a point in the body. We have thus achieved our goal that $\sig$ now describes the load in the body independent of which direction we cut the body in. With this information, let's use that to formulate the equilibrium equations!

## Equilibrium
In the previous example, we took the equilibrium to find the definition of the Cauchy stress. In the process, we saw that the the area terms where dominating over the volume terms as we went to an infitesimal volume. But let's do it again, now considering that the stress may vary over our small volume. For simplicity, we consider a regular hexahedron.
![](/assets/StressCube.svg)

The traction $\tv{t}_1$ will change slightly as we move from $x_0$ to $x_1=x_0+\dif x$, but we consider it's average value in the $y$ and $z$ directions. Similarly for $\tv{t}_2$ and $\tv{t}_3$ for $x,z$ and $x,y$ respectively. Taking the equilibrium equations for a volume load $\tv{b}$ and neglecting dynamic forces (i.e. quasi-static conditions) in the direction $\tv{e}_i$ we obtain
\begin{align*}
\left[\left[(\tv{t}_1(x_0)+\tv{t}_1(x_1)\right] \dif y \dif z + \left[(\tv{t}_2(y_0)+\tv{t}_2(y_1)\right] \dif x \dif z + \left[(\tv{t}_3(z_0)+\tv{t}_1(z_1)\right] \dif x \dif y\right] \cdot \tv{e}_i = -\tv{b} \cdot \tv{e}_i \dif V
\end{align*}
Inserting Cauchy's Theorem, and noting that the normal vectors are $\pm \tv{e}_j$, we get
\begin{align*}
\left[
\tv{e}_1 \cdot \left[\sig(x_1)-\sig(x_0)\right] \dif y \dif z + 
\tv{e}_2 \cdot \left[\sig(y_1)-\sig(y_0)\right] \dif x \dif z + 
\tv{e}_3 \cdot \left[\sig(z_1)-\sig(z_0)\right] \dif x \dif y
\right] \cdot \tv{e}_i = -\tv{b} \cdot \tv{e}_i \dif V
\end{align*}
where $\dif V = \dif x \dif y \dif z$. Changing to full index notation we get
\begin{align*}
\left[\sigma_{1i}(x_1)-\sigma_{1i}(x_0)\right] \dif y \dif z + \left[\sigma_{2i}(y_1)-\sigma_{2i}(y_0)\right] \dif x \dif z + \left[\sigma_{3i}(z_1)-\sigma_{3i}(z_0)\right] \dif x \dif y = -b_i \dif V
\end{align*}
Dividing by $\dif V=\dif x \dif y \dif z$, we obtain
\begin{align*}
\frac{\sigma_{1i}(x_0+\dif x)-\sigma_{1i}(x_0)}{\dif x} + \frac{\sigma_{2i}(y_0+\dif y)-\sigma_{2i}(y_0)}{\dif y} + \frac{\sigma_{3i}(z_0+\dif z)-\sigma_{3i}(z_0)}{\dif z} = -b_i
\end{align*}
Letting $\dif x \rightarrow 0$, $\dif y \rightarrow 0$, and $\dif z \rightarrow 0$, we obtain
\begin{align*}
\diff{\sigma_{1i}}{x} + \diff{\sigma_{2i}}{y} + \diff{\sigma_{3i}}{z} = \diff{\sigma_{ji}}{x_j} = -b_i
\end{align*}
which we can identify as the divergence of $\sig\trans$:
\begin{align}
\div{\sig\trans} + \tv{b} = \tv{0}
\end{align}
And this is the force equilibrium equation for a continuum. 

## Symmetry
Finally, we can show that $\sig$ is a symmetric tensor. By using Cauchy's Theorem (Equation \eqref{eq:CauchysTheorem}), we can give the tractions directly from the stress components in the following figure

![](/assets/StressSquare.svg)

with the same coordinate system as before. 

We then check that the counterclockwise moment around the $\tv{e}_3$ axis at $\tv{x}_0$ is zero. We denote the areas of the sides $\dif A_x = \dif y \dif z$ and $\dif A_y = \dif x \dif z$ for clarity.
\begin{align}
0 &= \sigma_{12} \dif A_x \dif x  - \sigma_{21} \dif A_y \dif y  \\
 &+ (\sigma_{22}(y_1)-\sigma_{22}(y_0)) \dif A_y \frac{\dif y}{2} - (\sigma_{11}(x_1)-\sigma_{11}(x_0)) \dif A_x \frac{\dif x}{2} \\
 &+ b_2 \dif V \frac{\dif x}{2} - b_1 \dif V \frac{\dif y}{2} 
\end{align}
where $b_1$ and $b_2$ are the $x$ and $y$ components of the volume load. Dividing by $\dif V=\dif x \dif y \dif z = \dif A_x \dif x = \dif A_y \dif y$ we obtain
\begin{align}
0 = \sigma_{12} - \sigma_{21} + \frac{\sigma_{22}(y_1)-\sigma_{22}(y_0)}{2} - \frac{\sigma_{11}(x_1)-\sigma_{11}(x_0)}{2} + \frac{b_2 \dif x}{2} - \frac{b_1 \dif y}{2}
\end{align}
Letting the size go to zero, $\sigma_{22}(y_1) \rightarrow \sigma_{22}(y_0)$, $\sigma_{11}(x_1) \rightarrow \sigma_{11}(x_0)$, $\dif x \rightarrow 0$, and $\dif y \rightarrow 0$, and we have
\begin{align}
\sigma_{12} = \sigma_{21}
\end{align}

Doing the same for the $xz$ and $yz$ planes, give the equivalent results, and we see that $\sigma_{ij}=\sigma_{ji}$, i.e. that $\ts{\sigma}=\tst{\sigma}$ showing that $\sig$ is a symmetric second order tensor. 
