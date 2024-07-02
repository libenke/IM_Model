# IM_Model
IM (Ianniruberto, G.and Marrucci, G.) Model  
### The Classic Single Mode IM Differential Model, or DCR-CS Model
ref: Ianniruberto, G.; Marrucci, G. A simple constitutive equation for entangled polymers with chain stretch. *Journal of Rheology* __2001__, 45, 1305-1318.  
Define an effective time $\tau_{eff}$ to replace the disengagement time $\tau_d$ in DE model:  
The classical DE model of $\mathbf{S}(t)$ is replaced to be:  

$$\begin{align}
&\overset{\bigtriangledown }{\mathbf{S}^2} +2\mathbf{S}^2(\boldsymbol{\kappa}:\mathbf{S})+\frac{2}{\tau}\mathbf{S}\cdot\left(\mathbf{S}-\frac{1}{3}\mathbf{I}\right) = \mathbf{0}\\
&\tau = \frac{1}{2\left(\frac{1}{\tau_d}+\boldsymbol{\kappa}:\mathbf{S}\right)} + \tau_R\\
&\frac{d\lambda}{dt} = \lambda\boldsymbol{\kappa}:\mathbf{S} - \frac{\lambda F(\lambda)-1}{\tau_R} \\
&F(\lambda) = \left(\frac{\lambda_{max}^2-\frac{\lambda^2}{3}}{\lambda_{max}^2-\lambda^2}\right)\left(\frac{\lambda_{max}^2-1}{\lambda_{max}^2-\frac{1}{3}}\right) \\
&\boldsymbol{\sigma} = 3G_N^0 F(\lambda)\lambda^2\boldsymbol{S}
\end{align}$$  

Here:  
$$
\overset{\bigtriangledown }{\mathbf{S}^2} = \mathbf{S}\cdot \dot{\mathbf{S}} + \dot{\mathbf{S}}\cdot \mathbf{S} -\boldsymbol{\kappa}\cdot\mathbf{S}^2-\mathbf{S}^2\cdot\boldsymbol{\kappa}^T
$$

## The Classical Single Mode IM Intergral Model
ref: Ianniruberto, G.; Marrucci, G. A simple constitutive equation for entangled polymers with chain stretch. *Journal of Rheology* __2001__, 45, 1305-1318.  
Define an effective time $\tau$ to replace the disengagement time $\tau_d$ in DE model:  
The classical DE model of $\mathbf{S}(t)$ is replaced to be:  

$$\begin{align}
\mathbf{S}(t) &= \int_{-\infty}^t \left[\frac{1}{\tau(t')}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,t')\right] dt' \\
\tau &= \frac{1}{2\left(\frac{1}{\tau_d}+\boldsymbol{\kappa}:\mathbf{S}\right)} +\tau_R \\
\frac{d\lambda}{dt} &= \lambda\boldsymbol{\kappa}:\mathbf{S} - \frac{\lambda F(\lambda)-1}{\tau_R} \\
F(\lambda)&=\left(\frac{\lambda_{max}^2-\frac{\lambda^2}{3}}{\lambda_{max}^2-\lambda^2}\right)\left(\frac{\lambda_{max}^2-1}{\lambda_{max}^2-\frac{1}{3}}\right) \\
\boldsymbol{\sigma} &= 3G_N^0F(\lambda)\lambda^2\boldsymbol{S}
\end{align}$$ 

## The Classic Multi Mode IM Integral Model
ref: Costanzo, S.; Huang, Q.; Ianniruberto, G.; Marrucci, G.; Hassager, O.; Vlassopoulos, D. Shear and Extensional Rheology of Polystyrene Melts and Solutions with the Same Number of Entanglements. *Macromolecules* __2016__, 49, 3925-3935.  
$$\begin{align}
\mathbf{S}_i(t) &= \int_{-\infty}^t \left[\frac{1}{\tau_i(t')}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau_i(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,t')\right] dt' \\
\tau_i(t) &= \frac{1}{2\left(\frac{1}{\tau_{i,eq}}+\boldsymbol{\kappa}:\mathbf{S}_i\right)} + \tau_R\\
\frac{d\lambda}{dt} &= \lambda\boldsymbol{\kappa}:\overline{\mathbf{S}} - \frac{\lambda F(\lambda)-1}{\tau_R} \\
F(\lambda)&=\left(\frac{\lambda_{max}^2-\frac{\lambda^2}{3}}{\lambda_{max}^2-\lambda^2}\right)\left(\frac{\lambda_{max}^2-1}{\lambda_{max}^2-\frac{1}{3}}\right) \\
\overline{\mathbf{S}} &= \int_{-\infty}^t \left[\frac{1}{\tau_d(t')}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau_d(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,t')\right] dt'\\
\tau_d(t) &= \frac{\sum_iG_i\tau_i^2(t)}{\sum_i G_i\tau_i(t)}\\
\boldsymbol{\sigma} &= C_Q F(\lambda)\lambda^2\sum{G_i\boldsymbol{S}_i}
\end{align}$$  
Here,
$C_Q = 6$ if $\mathbf{Q} = \frac{\mathbf{B}^{1/2}}{\mathrm{Tr}\mathbf{B}^{1/2}}$

# IM Tumbling Multi Mode Model 
ref: Costanzo, S.; Huang, Q.; Ianniruberto, G.; Marrucci, G.; Hassager, O.; Vlassopoulos, D. Shear and Extensional Rheology of Polystyrene Melts and Solutions with the Same Number of Entanglements. *Macromolecules* __2016__, 49, 3925-3935.
$$\begin{align}
\mathbf{S}_i(t) &= \int_{-\infty}^t \left[\frac{1}{\tau_i(t')}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau_i(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,t')\right] dt' \\
\tau_i(t) &= \frac{1}{2\left(\frac{1}{\tau_{i,eq}}+\boldsymbol{\kappa}:\mathbf{S}_i\right)} + \tau_R\\
\frac{d\lambda}{dt} &= \lambda\boldsymbol{\kappa}:\overline{\mathbf{S}} - \frac{\lambda F(\lambda)-1}{\tau_R} \\
F(\lambda)&=\left(\frac{\lambda_{max}^2-\frac{\lambda^2}{3}}{\lambda_{max}^2-\lambda^2}\right)\left(\frac{\lambda_{max}^2-1}{\lambda_{max}^2-\frac{1}{3}}\right) \\
\overline{\mathbf{S}} &= \int_{-\infty}^t \left[\frac{1}{\tau_d(t')}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau_d(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,t')\right] dt'\\
\tau_d(t) &= \frac{\sum_iG_i\tau_i^2(t)}{\sum_i G_i\tau_i(t)}\\
\boldsymbol{\sigma} &= C_Q F(\lambda)\lambda^2\sum{G_i\boldsymbol{S}_i}
\end{align}$$  
Here,
$C_Q = 6$ if $\mathbf{Q} = \frac{\mathbf{B}^{1/2}}{\mathrm{Tr}\mathbf{B}^{1/2}}$  
Tumbling term of $\lambda$ (shear case only):  
$$\begin{align}
\phi (t) &= \cos(2\pi \omega t)\exp (-\beta t)\\
\frac{d\lambda}{dt} &= \boldsymbol{\kappa}:\overline{\mathbf{S}}\phi(t)\lambda - \frac{F(\lambda)\lambda -1}{\tau_R}\\
\omega &= \frac{Wi_R^{-0.2}}{8\pi}\dot\gamma \\
\beta &= \frac{Wi_R^{-0.2}}{8}\dot\gamma
\end{align}$$

##### For step shear:
$$\begin{align}
\mathbf{S}(t) =& \int_{0}^t \left[\frac{1}{\tau(t')}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,t')\right] dt' \\
&+\int_{-\infty}^0 \left[\frac{1}{\tau(0)}\right] \exp \left[-\int_{t'}^t \frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,0)\right] dt' 
\end{align}$$
The second parts becomes:  
$$\begin{align}
\mathbf{S}(t) =& ...\\
&+\int_{-\infty}^0 \frac{1}{\tau(0)} \exp \left[-\int_{t'}^{0} \frac{dt''}{\tau(0)}-\int_0^t\frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,0)\right] dt' \\
=& ...\\
&+\int_{-\infty}^0 \frac{1}{\tau(0)} \exp \left[\frac{t'}{\tau(0)}-\int_0^t\frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,0)\right] dt' \\
=& ...\\
&+ \exp \left[\frac{t'}{\tau(0)}-\int_0^t\frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,0)\right] \bigg |_{t'=-\infty}^{t'=0}\\
=& ...\\
&+ \exp \left[-\int_0^t\frac{dt''}{\tau(t'')}\right] \mathbf{Q}\left[\mathbf{E}(t,0)\right]
\end{align}$$

#### Approximation of nonlinear strain measure $\mathbf{Q}\left[\mathbf{E}(t,t')\right]$
1. A highly accurate approximation for this strain measure is:  
ref: chapter 11.3 of Dealy, J. M. Larson, R. G. Read, D. J. "Structure and rheology of molten polymers, from structure to flow behavior and back again." __2018__, Carl Hanser Verlag GmbH Co KG.
$$\mathbf{Q} \approx  \left(\frac{5}{J-1}\right)\mathbf{B} - \left[\frac{5}{(J-1)(I_2+13/4)^{1/2}}\right]\mathbf{C}$$
here:  
$$J\equiv I_1+2(I_2+13/4)^{1/2}$$
$$I_1\equiv \mathrm{Tr}(\mathbf{B})$$
$$I_2\equiv \mathrm{Tr}(\mathbf{C})$$
$\mathbf{B}$ is the Finger tensor, and $\mathbf{C}$ is the Cauchy tensor.  
2. Marrucci et al. proposed another approximation of this strain measure, which gives a much improved prediction of the normal stress ratio in shear, namely $-N_2/N_1=1/4$ in the limit of small strains, as compared to Doi-Edwards value of 1/7. ref: Milner, S. T. Improved model of nonaffine strain measure. *Journal of Rheology* __2001__, 45, 1023-1028.
$$
\mathbf{Q}(\mathbf{E})=\frac{\mathbf{C}^{-1/2}}{\mathrm{Tr}(\mathbf{C^{-1/2}})}
$$
here:
$\mathbf{C}^{-1}=\mathbf{E}\cdot\mathbf{E}^T$ is the finger tensor. $\mathbf{E}$ is the inverse of displacement gradient tensor.

#### Finger Tensor $\mathbf{B}$ and Cauchy Tensor $\mathbf{C}$
For simple shear:  
$$\mathbf{F} (t_0,t_1)=\begin{pmatrix}
 1 & \gamma (t_1)-\gamma (t_0) & 0\\
 0 & 1 & 0\\
 0 & 0 &1
\end{pmatrix}$$
$$\mathbf{E} (t_0,t_1)=\begin{pmatrix}
 1 & \gamma (t_0)-\gamma (t_1) & 0\\
 0 & 1 & 0\\
 0 & 0 &1
\end{pmatrix}$$  
The Cauchy Tensor, also known as right Cauchy–Green deformation tensor, $\mathbf{C}$:  
$$
\begin{align}
\mathbf{C} &= \mathbf{F}^T\cdot\mathbf{F} \\
&=\begin{pmatrix}
 1 & \gamma & 0\\
 \gamma & 1+\gamma^2 & 0\\
 0 & 0 &1
\end{pmatrix}
\end{align}
$$
The Finger Tensor, also known as left Cauchy–Green deformation tensor, $\mathbf{B}$:  
$$\begin{align}
\mathbf{B} &= \mathbf{F}\cdot\mathbf{F}^{T} \\
&=\begin{pmatrix}
 1+\gamma^2 & \gamma & 0\\
 \gamma & 1 & 0\\
 0 & 0 &1
\end{pmatrix}
\end{align}$$
Sometimes, Finger strain tensor is also defined by inverse of Cauchy tensor as (ref: Milner, S. T. *Journal of Rheology* __2001__, 45, 1023-1028.):  
$$
\mathbf{B} = \mathbf{C}^{-1} = \mathbf{E}\cdot\mathbf{E}^T = \begin{pmatrix}
 1+\gamma^2 & -\gamma & 0\\
 -\gamma & 1 & 0\\
 0 & 0 &1
\end{pmatrix}
$$
Here, the diagonal is the same, while the $B_{0,1}, B_{0,2}, B_{1,2}$ parts are negative.  
