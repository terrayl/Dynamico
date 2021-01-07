# Dynamico
A set of dynamical adjustment codes

These codes perform a dynamical adjustment procedure as originally outlined in Deser et al. (2016). The main objective of the dynamical adjustment is to derive an estimate of the component of any physical variable variability due solely to atmospheric circulation changes. It has proven to be a very interesting approach in different attribution studies (Guo et al. 2019; Lehner et al. 2017; Lehner et al. 2018; Merrifield et al. 2017; O’Reilly et al. 2017). The method is based on (re-)constructed analogues of any atmospheric circulation variable (sea level pressure, geopotential height, atmospheric winds). It can be applied at the monthly, weekly or daily time scale. The original application (Deser et al. 2016) was based on monthly means. The original code has now also been adapted to daily time scale (Terray 2021). It can be used on both observational and model data. Note however that there are separate versions for model and observations as the basic ingredients and pre-processing phase do differ depending on the dataset nature (observed or simulated). The reader is referred to Deser et al. (2016) and Terray (2021) for further details.

The current code version is entirely written in NCL (NCAR Command Language, version 6.5.0, https://www.ncl.ucar.edu/). The codes are heavily commented to help users to easily follow the various methodological steps. I hope (fingers crossed) to be able to finalize a full python version later this year (2021!) when my IPCC work comes to an end ...

If you want to use any of the codes, please read the guidance file first where one can find some simple tips on how to proceed.

References:

Deser, C., Terray, L., and A. S. Phillips: Forced and internal components of winter air temperature trends over North America during the past 50 years: Mechanisms and implications. J. Climate, 29, 2237–2258, doi:10.1175/JCLI-D-15-0304.1, 2016

Guo, R., Deser, C., Terray, L., and Lehner, F.: Human influence on winter precipitation trends (1921–2015) over North America and Eurasia revealed by dynamical adjustment. Geophysical Research Letters, 46. https://doi.org/10.1029/2018GL081316, 2019

Lehner, F., Deser, C., Simpson, I., and Terray, L.: Attributing the US Southwest’s recent shift into drier conditions. Geophys. Res. Lett., 45. https://doi.org/10.1029/2018GL078312, 2018

Lehner, F., Deser, C., and Terray, L.: Toward a New Estimate of “Time of Emergence” of Anthropogenic Warming: Insights from Dynamical Adjustment and a Large Initial-Condition Model Ensemble. J. Climate, 30, 7739–7756, https://doi.org/10.1175/JCLI-D-16-0792.1, 2017

Merrifield, A., Lehner, F., Xie, S.-P., and Deser, C.: Removing circulation effects to assess central U.S. land-atmosphere interactions in the CESM Large Ensemble. Geophys. Res. Lett., 44, 9938–9946. https://doi.org/10.1002/2017GL074831, 2017

O’Reilly, C. H., Woollings, T., and Zanna, L.: The dynamical influence of the Atlantic multidecadal oscillation on continental climate. J. Climate, 30, 7213–7230, https://doi.org/10.1175/JCLI-D-16-0345.1, 2017

Terray, L.: A dynamical adjustment perspective on extreme event attribution. Weather and Climate Dynamics, to be submitted, 2021
