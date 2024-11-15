# CFCFluxes

Calculate air-sea fluxes of CFC-11 and CFC-12 according to OMIP protocols (Orr et al., 2017) and Dutay (2002).

The CFC flux ($F$) is given by:

$F_{CFC} = k_w (CFC_{sat} - CFC_{ocn})$

where:
- $k_w$ is the air-sea gas transfer velocity
- $CFC_{sat}$ is the saturation concentration in equilibrium with the water-vapor-saturated atmosphere
- $CFC_{ocn}$ is the sea surface tracer concentration

### Gas Transfer Coefficient

The instantaneous gas transfer coefficient ($k_w$) follows Wanninkhof (1992):

$k_w = a \left(\frac{Sc}{660}\right)^{-1/2} u^2 (1 - f)$

where:
- $a = 0.251$ (fitted constant)
- $Sc = A + BT_c + CT_c^2 + DT_c^3 + ET_c^4$ (Schmidt number polynomial)
  - $T_c$ is surface temperature in °C
- $u$ is 10-meter wind speed
- $f$ is sea ice coverage fraction

### Saturation Concentration

For small variations in surface pressure, $CFC_{sat}$ can be approximated as:

$CFC_{sat} \approx \frac{P_a}{P^0_{CFC}} F x_{CFC}$

where:
- $F$ is the CFC solubility function (Warner and Weiss, 1985)
- $x_{CFC}$ is atmospheric CFC concentration in ppt
- $P_a$ is sea surface pressure
- $P^0_{CFC} = 1 \text{ atm}$ (reference pressure)

## References

- Dutay et al. (2002)
- Wanninkhof (1992)
- Warner and Weiss (1985)
- Orr et al. (2017)
