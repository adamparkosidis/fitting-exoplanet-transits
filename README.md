#  Fitting TESS Exoplanet Transits with Python’s Numerical Libraries 

## Introduction

The two text files associated with this project contain time-series data for two stars with transiting exoplanets, observed by NASA’s Transiting Exoplanet Survey Satellite (TESS), which was launched in April 2018. During its initial 2-year survey mission, TESS has used its array of wide-field cameras to monitor 200.000 of the brightest stars near the Earth to look for small dips in the light from the star that correspond to (exo)planets transiting the face of the star and blocking some of the light. The two ‘light curves’ of TESS target stars provided here are given as time series with the first column showing the time in days (the values are given as Barycentred Julian Date – 2457000, although this is not relevant for this project) and the second column showing the calibrated photoelectron count rates for the star (as electrons/s although again, the exact flux units are not relevant for this Assignment). The light curves each cover about a month and the individual data points correspond to 2-minute intervals. There are gaps of a few days in the middle of the observations due to observing constraints on the satellite.

## Code overview

1. We read in the data using ```genfromtxt``` and ‘clean’ the data by removing rows from the resulting array which contain nan values (note that both the time and flux columns may contain these values). We plot the cleaned light curve.

2. We use the cleaned data to do the following analyses to  confirm the period of the transiting exoplanet:
    a. ```Using astropy.timeseries.LombScargle```, we measure and plot a Lomb-Scargle Periodogram of the time-series, to show at what potential
    periods the light curve shows the strongest variations.
    b. We smooth the light curve by averaging each data point together with the two data points on either side (i.e. we measure the average of ‘sliding       window’ of five data points centered on each data point). We plot the smoothed light curve along with ‘zoomed-in’ segment of the light curve showing     an eclipse. We use the smoothed light curve to confirm which peak is the likely orbital period of the transiting exoplanet.

3. Finally, we fit a model to the shape of the eclipse in the light curve, to constrain its parameters and measure a more accurate orbital period:
    a. We write a function which will calculate the following approximate functional form of an eclipse, when given the x-axis values:

    ![Eclipse shape](/images/eclipse_shape.png)

    Where `Δt_{total}` is the total duration of the eclipse (from start to finish), `Δt_{in}`  is the ingress time, i.e. the time from the start of the       eclipse until the maximum depth (`d`) of the eclipse is reached, `t_{cent}` is the time corresponding to the mid-point of the eclipse and `f_{star}`     is the flux level of the star when not eclipsed. You assume that the eclipse egress time (time of rising flux at the end of the eclipse) is the same     as the ingress time.

    b. We use the function with a subset of the cleaned (but unsmoothed) light curve, to fit one of the eclipses using ```scipy.optimize.curve_fit``` and     obtain the fit parameters described above. We plot the data and the model together on the same plot.
