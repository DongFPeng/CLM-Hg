# CLM-Hg

CLM-Hg is a modified version of the Community Land Model (CLM5) developed for representing terrestrial mercury cycling and mercury delivery from land to river systems.

This repository is maintained as a component of the **MOSART-Hg-WM** modeling framework based on CESM2.1.3. It contains the terrestrial mercury processes and model modifications used to provide mercury inputs to the coupled riverine mercury model.

## Important Notice

**This repository is provided solely for peer review and scientific evaluation.**

The source code, model modifications, parameterizations, configuration information, documentation, unpublished implementation details, figures, and other materials contained in this repository are not released for unrestricted public reuse.

Any reuse, redistribution, adaptation, incorporation into other models or software, derivative development, or use of unpublished materials from this repository requires **prior written permission** from **Dong Peng or the corresponding author(s)** of the associated work.

If you wish to reuse or reference repository-specific materials, model implementation details, parameterizations, or unpublished methods, please contact the authors before use.

These restrictions apply to the authors' original modifications and materials and do not supersede the licenses or terms governing CESM, CLM, or other third-party software.

## Model base

```text
CESM version: 2.1.3
Land model:   CLM5
Component:    CLM-Hg
```

CLM-Hg extends the terrestrial component of CESM to represent processes relevant to the mobilization and transfer of mercury from land to river networks.

## Main developments

The model framework includes modifications related to:

- terrestrial mercury storage and redistribution
- atmospheric mercury deposition to land
- mercury mobilization associated with runoff
- soil erosion and particulate mercury transport
- dissolved and particulate mercury delivery to rivers
- anthropogenic mercury inputs
- legacy mercury stored in terrestrial reservoirs
- source attribution of land-to-river mercury fluxes
- climate- and land-use-driven changes in mercury mobilization

The model is designed to provide physically consistent terrestrial mercury inputs for subsequent routing and transformation within MOSART-Hg.

## Role in MOSART-Hg-WM

Within the coupled framework:

```text
Atmosphere / Anthropogenic inputs
              ↓
           CLM-Hg
              ↓
     Land-to-river Hg fluxes
              ↓
          MOSART-Hg
              ↓
     River transport / retention
              ↓
             Ocean
```

CLM-Hg provides terrestrial mercury fluxes associated with runoff, erosion, and other land-surface processes to the river network.

The coupled model configuration is maintained in:

**MOSART-Hg-WM**  
`https://github.com/DongFPeng/MOSART-Hg-WM`

## Repository use

This repository is primarily intended to preserve the specific CLM-Hg development history and model versions used in MOSART-Hg-WM.

For coupled simulations, users should refer to the parent MOSART-Hg-WM repository, which records the exact CLM-Hg and MOSART-Hg commits used for each model configuration.

## Publications

The development and application of the coupled terrestrial-river mercury framework are described in the following studies:

1. **Peng, D. et al. (2025). Human perturbations to mercury in global rivers. _Science Advances_, 11, eadw0471.**  
   DOI: `10.1126/sciadv.adw0471`

2. **Peng, D. et al. (2026). Tracing Mercury from Land to River: Global Sources, Retention, and Implications for Sustainability. _Environmental Science & Technology_, 60, 21885–21900.**  
   DOI: `10.1021/acs.est.6c03367`

Additional publications and research outputs are available through the author's ResearchGate profile.

## Data and model input

Large CESM input datasets, atmospheric forcing data, NetCDF input files, restart files, model outputs, and temporary processing files are not maintained in this repository.

Required CESM input datasets should be obtained separately and configured according to the local computing environment.

## License and use restrictions

The underlying CESM and CLM source code remains subject to its original license and terms of use.

The authors' original mercury-related modifications, parameterizations, configurations, documentation, and unpublished implementation details are provided for **peer-review and scientific-evaluation purposes only**, unless explicit written permission for further use has been granted.

Permission requests should be directed to Dong Peng or the corresponding author(s) of the associated publication.

## Contact

**Dong Peng**  
Email: dongpeng506@gmail.com  
ResearchGate: `https://www.researchgate.net/profile/Dong-Peng-4`

For questions regarding CLM-Hg, requests for reuse, or permission to use repository-specific materials, please contact the authors before use.


# CTSM


The Community Terrestrial Systems Model.

This includes the Community Land Model (CLM5.0 and CLM4.5) of the Community Earth System Model.

For documentation, quick start, diagnostics, model output and
references, see

http://www.cesm.ucar.edu/models/cesm2.0/land/

and

https://escomp.github.io/ctsm-docs/

For help with how to work with CTSM in git, see

https://github.com/ESCOMP/ctsm/wiki/Getting-started-with-CTSM-in-git

and

https://github.com/ESCOMP/ctsm/wiki/Recommended-git-setup

To get updates on CTSM tags and important notes on CTSM developments
join our low traffic email list:

https://groups.google.com/a/ucar.edu/forum/#!forum/ctsm-dev

(Send email to ctsm-software@ucar.edu if you have problems with any of this)
