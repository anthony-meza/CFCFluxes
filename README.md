# CFCSolubility

The purpose of this code is to faciliate the easily 
calculation of air-sea fluxes of CFC-11 and CFC-12. The flux, 
$Q$ is given by, 
$$Q = k \left(F P_{cfc} \frac{P}{P_0} - C_s \right)$$, 
according to Dutay (2002). Here, 
- $k$ is the air-sea gas transfer velocity
- $F$ is the CFC solubility function from Warner and Weiss (1985)
- $P_{CFC}$ atmospheric partial pressure in dry air at one atmosphere total pressure. It has units atm. 
- $P$ is the sea surface pressure. It has units in atm
- $P_0 = 1 \text{ atm}$
- $C_s$ is the sea surface tracer concentration 

Typically, CFC data is given in atmospheric concentration, $C_a$. If $C_a$ has units parts per trillion, then $P_{cfc} \approx = C_a \times 10^{12} / 1 \text{ atm}$.
