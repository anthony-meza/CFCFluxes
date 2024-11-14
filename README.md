# CFCSolubility

The purpose of this code is to faciliate the easily 
calculation of air-sea fluxes of CFC-11 and CFC-12. The flux, 
$Q$ defined in the OMIP protocols and is given by, 

$F_{CFC} = k_w (CFC_{sat} - CFC_{ocn})$  
according to Dutay (2002). Here, 
- $k_w$ is the air-sea gas transfer velocity
- $CFC_{sat}$ is the saturation concentration in equilibrium
with the water-vapor-saturated atmosphere at a total atmospheric pressure Pa
- $CFC_{ocn}$ is the sea surface tracer concentration 

The instantaneous gas tranfer coefficient $k_w$ is parameterization from Wanninkhof (1992) 

$$k_w = a \left(\frac{Sc}{660}^{-1/2}\right) u^2 (1 - f).$$
Here, 
- $a = 0.251$ is a fitted constant \
- $Sc = A + B T_c + C T_c^2 + D T_c^4 + E T_c^4$ is the polynomial Schmidt number where $T$ is the surface temperature in deg Celsius.]
- $u$ is the 10-meter wind speed 
- $f$ is the fraction of sea ice coverage  

Assuming that surface pressure only varies slightly, $CFC_{sat}$ may be 
approximated as: 
$$ CFC_{sat} \approx = \frac{P_{CFC}}{P_{CFC}^0} F x_{CFC}$$

- $F$ is the CFC solubility function from Warner and Weiss (1985)
- $x_{CFC}$ are atmospheric records of CFC-11 and CFC-12 in parts per trillion (ppt)
- $P_a$ is the sea surface pressure
- $P^0_{CFC} = 1 \text{ atm}$ is the reference 
