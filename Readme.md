# Boussinesq rotating thermal convection(RTC) in a full sphere

The model equations are

$$
(\partial_t -\Delta){\bf u}= {\bf u}\times(\nabla\times{\bf u})- \frac{1}{E}\hat{\bf z}\times{\bf u}  -\nabla \Pi + \frac{Ra}{E} \Theta{\bf r}\\
$$

$$
(\partial_t -\frac{1}{Pr}\Delta)\Theta = -{\bf u}\cdot\nabla\Theta + {\bf u}\cdot{\bf r}
$$

with the parameters defined as

$$
Pr = \frac{\nu}{\kappa};\quad E=\frac{\nu}{2\Omega r_o^2};\quad Ra=\frac{\gamma \alpha \beta r_o^4}{2\nu\Omega}
$$

where $\nu$, $\kappa$, $\alpha$ are the kinematic viscosity, thermal conductivity, and thermal expansion coefficents; $\Omega$ is the rotation rate and $r_o$ is the radius of the sphere; $\gamma$ and $\beta$ are defined via the gravitational acceleration, ${\bf g} = -\gamma {\bf r}$, and the background temperature profile for internal heating, $\nabla\Theta_b = -\beta {\bf r}$. The equations have been non-dimensionalised via the viscous time-scale, and using $\beta r_o^2$ as temperature scale.


The system is solved using a Toroidal/Poloidal decomposition of the velocity

$$
{\bf u} = \nabla\times T {\bf r} + \nabla\times\nabla\times P {\bf r} 
$$
