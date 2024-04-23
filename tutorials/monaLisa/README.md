## Mona Lisa Tutorial
This case replicates the Mona Lisa spot melt pattern reported by Plotkowski et al. [1]. The case has the properties below.

### Material
- IN625 used in tutorial (IN718 reported in [1])

### Beam Properties
- Power: 720W (12 mA melt current, 60 kV electron source)
- Size: 200um beam width (100 um D2sigma, estimated 30um depth)
- Absorption: 90% (much higher for EB-PBF than L-PBF)
- Used standard superGaussian heat source model

### Scan Path
The scan path was adopted from the scan path A reported by [1]

### Mesh
The scan path covers an area of approximately 35 mm x 55 mm. The mesh has been extended 1 mm to either side of the scan path. A base mesh size of 80 um is used, and adaptive mesh refinement with 3 refinement levels is used to achieve a 10 um cell size in the scan path. The resource-optimized algorithm is used to target a mesh size of 10,000 cells per core.


### References
[1] Plotkowski, A., Ferguson, J., Stump, B., Halsey, W., Paquit, V., Joslin, C., Babu, S. S., Marquez Rossy, A., Kirka, M. M., & Dehoff, R. R. (2021). A stochastic scan strategy for grain structure control in complex geometries using electron beam powder bed fusion. Additive Manufacturing, 46. https://doi.org/10.1016/j.addma.2021.102092
