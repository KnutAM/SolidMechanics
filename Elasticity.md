+++
 title = "Elasticity"
+++

**Intended learning outcomes**
* Define what an elastic material model is
* Name a other types of material behavior
* Define the linear elastic material model and the hyperelastic material model
* Show if a linear elastic material law is isotropic or not
* Know the properties of key classes of anisotropy: Transverse isotropy and orthotropy


# Elasticity
The definition of an elastic material model is that the stress, $\sig$, is a unique function of the strain, $\eps$. That is, if we know $\eps$, we can calculate $\sig$ via the function $\sig(\eps)$. 
Consider the following figures that show different material behaviors
![](/assets/MaterialBehavior.svg)

~~~
<ol type="a">
  <li>The stress-strain relationship is nonlinear, but for every strain value, there is only one stress value: the material behavior shown is elastic.</li>
  <li>The stress-strain relationship is nonlinear, including softening. The same stress value occurs for two different strains, but for each strain value, only a single stress value exists: the material behavior shown is elastic.</li>
  <li>After loading the material, the stress-strain curve follows a different path back. The behavior seen here is typical for a damaged material, when the stiffness is reduced due to microscopic cracks in the material. Such behavior is not purely elastic. </li>
  <li>Also here the stress-strain curve follows a different path upon unloading. Here, the tangent is initially the same as the first part of the loading curve. Such behavior is typical for plastic material models, and the shown behavior is not purely elastic. </li>
</ol>
~~~

Another class of material behavior not shown above, are rate dependent materials. Typically, polymers have a viscoelastic behavior, and for the same strain we get a higher stress if we increase the strain rate (i.e. the loading rate). If a material is rate dependent, we cannot determine $\sig$ from $\eps$ alone (we also need $\dot{\eps}$), hence this is not an elastic material model. 

In general, many materials behave (approximately) elastic within some limits. For example, metals are elastic when the stress is below the yield limit and polymers behave almost elastic if the loading rate is sufficiently slow and the strain is sufficently low. 

## Hyperelasticity
Hyperelasticity is defined by the existence of a so-called potential $\varPsi(\eps)$ (a scalar function of the strain), such that 
\begin{align}
\sig = \pdiff{\varPsi}{\eps} \label{eq:hyperelasticity}
\end{align}



## Linear elasticity
The most general form of linear elasticity is given by 
\begin{align}
\sig = \tf{E}:\eps \label{eq:linearelasticity}
\end{align}
where $\tf{E}$ is the 4th order elastic stiffness tensor. This model is hyperelastic provided that $\tf{E}$ is major symmetric, with $\varPsi = 0.5\eps:\tf{E}:\eps$. This is seen by taking the derivative
\begin{align*}
\pdiff{\sig}{\eps} = \tf{E}
\end{align*}
from Equation \eqref{eq:linearelasticity}. Using Equation \eqref{eq:hyperelasticity}, we get that
\begin{align*}
\pdiff{\sig}{\eps} = \frac{\partial^2 \varPsi}{\partial \eps \partial \eps} = \tf{E}
\end{align*}
Considering this expression in index notation, we have
\begin{align*}
\tfind{E}{ijkl} = \frac{\partial^2 \varPsi}{\partial \epsilon_{ij} \partial \epsilon_{kl}} = \frac{\partial}{\partial \epsilon_{ij}}\left[\pdiff{\varPsi}{\epsilon_{kl}}\right] = \frac{\partial}{\partial \epsilon_{kl}}\left[\pdiff{\varPsi}{\epsilon_{ij}}\right]
\end{align*}
because we can choose which value to differentiate with first ($\epsilon_{ij}$ or $\epsilon_{kl}$). Therefore, we have that $\tfind{E}{ijkl}=\tfind{E}{klij}$, i.e. that $\tf{E}$ is major symmetric. Furthermore, since $\eps$ is symmetric, $\tf{E}$ must also be minor symmetric. In 3 dimensions, a 4th order tensor has $3^4=81$ components. But due to the minor symmetries of $\tf{E}$, this is reduced to $6^2=36$. Furthermore, with the major symmetry there are only 21 remaining independent components, which is the maximum number of possible elastic parameters for a linear elastic material. 

### Isotropic linear elasticity
The behavior of an isotropic material is independent of the loading direction. So given our function $\sig = \ts{f}(\eps)$ for an elastic material, an isotropic elastic material fulfills
\begin{align}
\ts{R}\sig(\eps)\tst{R} = \ts{f}\left(\ts{R}\eps\tst{R}\right)
\end{align}
where $\ts{R}$ is a proper orthogonal rotation tensor. We use the function $\ts{f}$ just to highlight the different between the value of $\sig$ and the function $\sig(\eps)$. 
In the case of linear elasticity, this relationship translates into
\begin{align}
\tf{E} = \left[\tst{R}\opu\tst{R}\right] : \tf{E} : \left[\ts{R}\opu\ts{R}\right]
\end{align}

Isotropic linear elasticity is given by
\begin{align}
\sig = 2G \eps\dev + 3K \eps\sph \label{eq:isotropic}
\end{align}
where $G$ is the shear modulus and $K$ is the bulk modulus. Hence, only 2 of the possible 21 elastic parameters are required. Applying the definition of isotropic material behavior directly, we obtain in index notation
\begin{align*}
R_{ij}\left[2G \epsilon_{jk} + \frac{3K-2G}{3}\epsilon_{mm} \delta_{jk}\right]R_{kl}\trans &= 2G R_{ij}\epsilon_{jk}R_{kl}\trans  +  \frac{3K-2G}{3} \left[\left[ R_{mn}\epsilon_{no}R_{om}\trans\right] \right] \delta_{il} \\
2G R_{ij}\epsilon_{jk}R_{kl}\trans + \frac{3K-2G}{3}\epsilon_{mm} R_{ij}\delta_{jk}R_{kl}\trans &= 2G R_{ij}\epsilon_{jk}R_{kl}\trans  +  \frac{3K-2G}{3} \left[\left[R_{om}\trans R_{mn}\epsilon_{no}\right] \right] \delta_{il} \\
2G R_{ij}\epsilon_{jk}R_{kl}\trans + \frac{3K-2G}{3}\epsilon_{mm} \delta_{il} &= 2G R_{ij}\epsilon_{jk}R_{kl}\trans  +  \frac{3K-2G}{3} \epsilon_{nn}  \delta_{il} \\
\end{align*}
where we use that $\ts{R}\tst{R}=\ts{I}$. This equality shows that Equation \eqref{eq:isotropic} describes a material behavior that is independent of the loading direction. 

Differentiating the linear elastic relation to obtain $\tf{E}$, we get
\begin{align}
\pdiff{\sig}{\eps} = 2G \left[\tf{I} - \frac{1}{3}\ts{I}\otimes\ts{I}\right] + K \ts{I}\otimes\ts{I}
\end{align}

### Linear transverse isotropic elasticity
Transverse isotropy implies that a material behaves isotropic for in-plane loading in a specific plane. An example of such a material, is a uniaxially fiber-reinforced material as illustrated below. Here, the behavior transverse to the fibers is isotropic. However, the material will be stiffer when loaded along the fibers compared to transverse to the fibers. 

![](/assets/TransverselyIsotropic.svg)

For the material symmetry shown above, this implies that the elastic relation is unaffected by rotations about the $\tv{n}_3$. More precisely
\begin{align}
\ts{R}_3\sig(\eps)\tst{R}_3 = \ts{f}\left(\ts{R}_3\eps\tst{R}_3\right), \quad \tst{R}_3\ts{R}_3 = \ts{I}, \; \ts{R}_3\tv{n}_3 = \tv{n}_3
\end{align}

For this specific symmetry, there are only 5 indpendent material parameters. 

### Linear orthotropic elasticity
While it is possible to have models for anisotropy with 21 independent material paramters as shown above, very few models consider this. 
Orthotropic material behavior is often the most general form of elastic anisotropy that is assumed. Formally, it is defined by invariance under so-called orthotropic transformations
\begin{align}
\ts{\hat{R}}\sig(\eps)\tst{\hat{R}} = \ts{f}\left(\ts{\hat{R}}\eps\tst{\hat{R}}\right), \quad \ts{\hat{R}} = \ts{I} - 2\tv{n}^{(i)}\otimes\tv{n}^{(i)}, \quad \tv{n}^{(i)}\cdot\tv{n}^{(j)}=\delta_{ij}
\end{align}
Here, $\tv{n}^{(i)}$ represents the three symmetry planes. Note, that the above relation must only hold for a specific set of symmetry planes for the material to be orthotropic. The transformation tensors fullfill the identity $\ts{\hat{R}}=\tsi{\hat{R}}$ and in a coordinate system aligned with $\tv{n}^{(i)}$, the possible $\ts{\hat{R}}$ are represented as
\begin{align*}
\begin{bmatrix}-1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1\end{bmatrix}, \quad \begin{bmatrix}1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 1\end{bmatrix} \quad \begin{bmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -1\end{bmatrix}
\end{align*}

But these definitions are hard to grasp the physical implications of. To better understand this, we consider two more physical explanations

1. There exists three special perpendicular directions. If we apply a normal load in any of these directions, no shear strains will arise. 
2. There exists three special perpendicular directions. Consider a shear strain on a plane perpendicular to one of these directions. If we reverse the shear strain but maintain all other strain components constant, the corresponding shear stress component reverses, while the remaining stress components remain unchanged. 

For this specific symmetry, there are 9 independent material parameters. It can be very challenging to determine these if we don't know the orthotropic directions. However, if those directions are known, it is fairly straight-forward to identify the parameters. 