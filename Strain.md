+++
 title = "Strain"
+++

**Intended learning outcomes**

* Explain how the deformation gradient $\ts{F}$ represents the deformation of a volume element
* Show that the small strain tensor, $\ts{\epsilon}$, is a reasonable strain measure for small rotations


# Strain
The strain is a local measure of the shape-changes of a structure. As it is a local measure, its connection with the stress is not related to the structure's shape. Therefore, it is suitable to model the material behavior as the stress as a function of the strain: This material behavior can be characterized independently from the structure. For example, we can perform a laboratory test on a test sample to measure the material's stress-strain relationship. We can then apply this relationship to determine the local stress states when simulating an entire structure. But the question then is, how to define such a local strain measure? We have two key requirements

1. It should describe how an infitesimal volume element (think cube) changes shape
2. It should be independent of rigid body motion (translation and rotation)

In addition, for convenience, we would like to have a strain measure such that if it is zero there is no deformation
3. It should be zero if there are no deformations. 


## The deformation gradient and the small strain tensor
The deformation gradient, $\ts{F}$, measures the *gradient* of the *deformed coordinates*, $\tv{x}=\tv{X}+\tv{u}$ (the name *displacement gradient* is sometimes used for $\ts{F}-\ts{I}=\pdiffil{\tv{u}}{\tv{X}}$)
\begin{align}
\ts{F} = \gradX{\tv{x}} = \pdiff{\tv{x}}{\tv{X}} = \pdiff{\tv{u}}{\tv{X}} + \ts{I}
\end{align}
Here, $\tv{x}$ are the deformed coordinates and $\tv{X}$ the undeformed coordinates. That is, if we have the displacements $\tv{u}$, then $\tv{x}=\tv{X}+\tv{u}$. In 2 dimensions, the effect of the deformation gradient on an area, which is initially square, can be illustrated as

![](/assets/DeformationGradient.svg)

Based on this figure, it is clear that the deformation gradient describes how an infitesimal volume element changes shape. But it is not zero when the strain is zero - it is the identity tensor. Defining the strain as $\ts{F}-\ts{I}$ solves this problem. However, consider a pure rotation of the square:

![](/assets/DeformationGradientRotation.svg)

For the small rotation, $\beta=5^\mathrm{o}$, we see that we obtain a change in the deformation gradient's off-diagonal (shear) terms. However, these are of equal magnitude but with opposite signs. So in the symmetric part of the deformation gradient is negligibly affected by the small rotation. So for small rotations, we can thus define a strain measure as

\begin{align}
\ts{\epsilon} = \left[\gradX{\tv{x}}\right]\sym - \ts{I} = \left[\pdiff{\tv{x}}{\tv{X}}\right]\sym  - \ts{I} = \left[\pdiff{\tv{u}}{\tv{X}}\right]\sym
\end{align}
This strain measure is the small strain tensor, and we will use it throughout the linear continuum mechanics course. 

Considering the larger rotation, $\beta=30^\mathrm{o}$, even the diagonal (normal) terms are affected. Therefore, $\ts{\epsilon}$ is not suitable for large rotations. 

In later courses dealing with nonlinear continuum mechanics, so-called finite strain measures can be introduced. Once example is the Green-Lagrange strain, $\ts{E}$, defined as
\begin{align}
\ts{E} = 0.5\left[\tst{F}\ts{F} - \ts{I}\right]
\end{align}

A current coordinates after a pure rotation are given by $\tv{x} = \ts{R}\tv{X}$, where $\ts{R}$ is a proper orthogonal rotation tensor. For this case,  we obtain $\ts{F} = \ts{R}$ and our two strain measures become
\begin{align}
\ts{\epsilon} &= 0.5 \left[\ts{R} + \tst{R}\right] - \ts{I} \neq \ts{0} \\
\ts{E} &= 0.5\left[\tst{R}\ts{R} - \ts{I}\right] = \ts{0}
\end{align}

So why do we not use this strain tensor instead? It involves square terms of the deformation gradient, and is therefore not linear. When solving structural problems, the linearity of the problem makes it much easier to solve. If we consider the linearization of the green-lagrange strain tensor, we would like to see how it is affected by a small deformation, i.e. $\ts{F} = \ts{I} + \delta\ts{F}$, where $\ts{\delta F}\ll\ts{I}$ we have
\begin{align*}
\ts{E} &= 0.5\left[\left[\tst{I}+\tst{\delta F}\right] \left[\ts{I}+\ts{\delta F}\right] - \ts{I} \right] = 0.5\left[\ts{I} + \ts{\delta F} + \tst{\delta F} + \tst{\delta F}\ts{\delta F} - \ts{I}\right] \\
&\approx \ts{\delta F}\sym = \left[\ts{I} + \ts{\delta F}\right]\sym - \ts{I} = \ts{\epsilon}
\end{align*}
where the fact that $\ts{\delta F}\ll \ts{I}$ causes the term $\tst{\delta F}\ts{\delta F}$ to be negligible compared to $\ts{\delta F}$. 
