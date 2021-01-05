# Dynamico
A program to perform a dynamical adjustment procedure

This code performs a dynamical adjustment procedure as outlined in Deser et al. (2016). The main objective of the dynamical adjustment is to derive an estimate of the component of any physical variable variability due solely to atmospheric circulation changes. It has proven to be a very interesting approach in different attribution studies (Guo et al. 2019; Lehner et al. 2017; Lehner et al. 2018). The method is based on (re-)constructed analogues of any atmospheric circulation variable (sea level pressure, geopotential height, atmospheric winds). It can be applied at the monthly, weekly or daily time scale. The original application (Deser et al. 2016) was based on monthly means and the code has now also been extended to daily time scale (Terray 2021). It can be used on both observational and model data. Note however that there are separate versions for model and observations as the basic ingredients and pre-processing phase do differ depending on the dataset nature (observed or simulated). The reader is referred to Deser et al. (2016) for further details.

The current code version is entirely written in NCL (NCAR Command Language, version 6.5.0). I hope (fingers crossed) to be able to finalize a full python version later this year when my IPCC work comes to an end ...

References:

Deser, C., Terray, L., and A. S. Phillips: Forced and internal components of winter air temperature trends over North America during the past 50 years: Mechanisms and implications. J. Climate, 29, 2237–2258, doi:10.1175/JCLI-D-15-0304.1, 2016

Guo, R., Deser, C., Terray, L., and Lehner, F.: Human influence on winter precipitation trends (1921–2015) over North America and Eurasia revealed by dynamical adjustment. Geophysical Research Letters, 46. https://doi.org/10.1029/2018GL081316, 2019

Lehner, F., Deser, C., Simpson, I., and Terray, L.: Attributing the US Southwest’s recent shift into drier conditions. Geophys. Res. Lett., 45. https://doi.org/10.1029/2018GL078312, 2018

Lehner, F., Deser, C., and Terray, L.: Toward a New Estimate of “Time of Emergence” of Anthropogenic Warming: Insights from Dynamical Adjustment and a Large Initial-Condition Model Ensemble. J. Climate, 30, 7739–7756, https://doi.org/10.1175/JCLI-D-16-0792.1, 2017

Terray, L.: A dynamical adjustment perspective on extreme event attribution. Weather and Climate Dynamics, to be submitted, 2021
