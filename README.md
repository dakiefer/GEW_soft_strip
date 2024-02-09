# GEW soft strip [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10372576.svg)](https://doi.org/10.5281/zenodo.10372576)

**Compute guided elastic waves (GEWs) in a soft strip.** 

This code computes dispersion curves of an anisotropic waveguide with rectangular cross-section. It uses the spectral collocation method (SCM) and is an extension to "[GEW dispersion script](https://github.com/dakiefer/GEW_dispersion_script)".

Three cases are included: 
1. `strip_elastic_SCM.m`: linear-elastic (possibly anisotropic) strip.
2. `strip_elastic_SA_SCM.m`: as before but allows to choose the symmetry of waves.
3. `strip_soft_viscoAE_SCM.m`: includes **acoustoelasticity** and **viscosity** for a soft, nearly-incompressible strip.

Acoustoelasticity and viscosity are inherently present in common soft matter, such as biological tissues. Acoustoelasticity describes the effect of pre-stress on the wave propagation, which we account for using a compressible Mooney-Rivlin hyperelastic material model. Viscoelastic losses are included with a fractional Kelvin-Voigt model. The difficulty relies in the interdependence of both effects. 

The methods have been presented in:

> A. Delory, D. A. Kiefer, M. Lanoy, A. Eddi, C. Prada, and F. Lemoult, “Viscoelastic dynamics of a soft strip subject to a large deformation,” *Soft Matter*, Jan. 2024, doi: [10.1039/D3SM01485A](https://doi.org/10.1039/D3SM01485A).

Code repository: [<img src="https://www.svgrepo.com/show/35001/github.svg" alt="GitHub" width="27px" />](https://github.com/dakiefer/gew_zgv_computation) [https://github.com/dakiefer/gew_soft_strip](https://github.com/dakiefer/gew_soft_strip)

## How to use

1. Change into the `GEW_soft_strip` folder or add it to the Matlab path.
2. Execute the desired script, i.e., `strip_elastic_SCM.m`, `strip_elastic_SA_SCM.m` or `strip_soft_viscoAE_SCM`. 
3. If you want to implement a different hyperviscoelastic material model, you can create it by adapting `initiateTensor.nb` and `createTensorForMatlab.nb` and running the latter in Mathematica. Copy the output of the last command as *plain text* and replace the marked block in `customTensor.m` with the code you copied. 

## Limitations 

Mixed boundary conditions (Dirichlet and Neumann) can lead to instabilities at the corners. Switch to another discretization method, e.g., Finite Elements, if this is needed.

## Mathematical background 

The SCM discretizes boundary-value problems based on the strong form. The two-dimensional implementation that we use is explained in 
> J. A. Weideman and S. C. Reddy, “A MATLAB Differentiation Matrix Suite,” ACM Trans. Math. Softw., vol. 26, no. 4, pp. 465–519, 2000, doi: [10.1145/365723.365727](http://doi.org/10.1145/365723.365727).

## Dependencies

`chebdif.m` from DMSUITE is bundled and can be found at

J.A.C Weideman (2022). DMSUITE (https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite), MATLAB Central File Exchange. Retrieved August 18, 2022.

`ToMatlab.m` is bundled from

Harri Ojanen (1999). Mathematica Expression to Matlab m-file Converter (https://library.wolfram.com/infocenter/MathSource/577/), Wolfram Library Archive. Retrieved  June 23rd, 2022.

## Citing this software

If this code is useful to you, please cite it as:

> D. A. Kiefer, A. Delory, and F. Lemoult. GEW soft strip (2023), doi [10.5281/zenodo.10372576](http://doi.org/10.5281/zenodo.10372576). https://github.com/dakiefer/GEW_soft_strip

together with the related publication:

> A. Delory, D. A. Kiefer, M. Lanoy, A. Eddi, C. Prada, and F. Lemoult, “Viscoelastic dynamics of a soft strip subject to a large deformation,” *Soft Matter*, Jan. 2024, doi: [10.1039/D3SM01485A](https://doi.org/10.1039/D3SM01485A).

## Authors

Code created 2022–2023 by

**Daniel A. Kiefer**, Institut Langevin, ESPCI Paris, Université PSL, France<br/>
[daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr) &nbsp; ● &nbsp; [dakiefer.net](https://dakiefer.net) &nbsp; ● &nbsp; [Google Scholar](https://scholar.google.de/citations?user=odSy3v4AAAAJ&hl=en) &nbsp; ● &nbsp; Follow me on [ResearchGate](https://www.researchgate.net/profile/Daniel-Kiefer-5)!

**Alexandre Delory**, Institut Langevin, ESPCI Paris, Université PSL, France<br/>
[delory.alexandre@free.fr](mailto:delory.alexandre@free.fr) &nbsp; ● &nbsp; [Google Scholar](https://scholar.google.de/citations?hl=en&user=OgjaLqIAAAAJ).

**Fabrice Lemoult**, Institut Langevin, ESPCI Paris, Université PSL, France<br/>
[fabrice.lemoult@espci.psl.eu](mailto:fabrice.lemoult@espci.psl.eu) &nbsp; ● &nbsp; [Google Scholar](https://scholar.google.de/citations?user=Gy6ImbgAAAAJ).

[<img src="https://user-images.githubusercontent.com/3725269/185571121-f5fcd518-32de-40b2-b4b1-f4ef0610ccd1.svg" alt="Logo Institut Langevin" width="250px" />](https://www.institut-langevin.espci.fr) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<img src="https://user-images.githubusercontent.com/3725269/185570398-ca2796ab-2bd3-4171-a7a6-af1f74014504.svg" alt="Logo ESPCI Paris" width="260px" />](https://www.espci.psl.eu/en/)
