import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
import os

#########################
### Define Functions ####
#########################

def data_cleaning(filename):
    data = np.genfromtxt(filename)
    data_cleaned=data[~np.isnan(data).any(axis=1)]

    fig = plt.figure()
    plt.plot(data[:,0],data[:,1])
    plt.xlabel('Time (BJD)')
    plt.ylabel('Flux (electorns/sec)')
    plt.title('Light Curve during Transit')
    plt.grid()
    return data_cleaned

def period_exoplanet(data_cleaned, which_txt):
    frequency_array = np.linspace(0.1,2,10000) # frequency values to evaluate the periodogram over
    power =LombScargle(data_cleaned[:,0],data_cleaned[:,1]).power(frequency_array) # the distribution of power into frequency components composing that signal

    '''
    The first subplot shows the power versus frequency, while the second one the power versus 1/frequency (Period).
    By seeing the second subplot we can see different peeks which correspond to different candidates for the period.
    Comparing each one of them with the light curve (smoothed data) we can conclude that the peek corresponding to 
    ~ 3.2 days is the period of the planet for the second data sheet. For the first data sheet, things are not that
    clear, but my guess is that the peek corresponding to ~ 1 day is the period of the planet.
    '''
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.tight_layout()
    ax1.plot(frequency_array,power)
    ax2.plot(1/frequency_array,power)
    ax1.grid()
    ax2.grid()
    
    
    kernel_size = 5
    kernel = np.ones(kernel_size) / kernel_size    # define a kernel to convolve the data
    flux_smoothed=np.convolve(data_cleaned[:,1],kernel, mode='same') # the smoothed data after the convolution
    
    '''
    In case that the "#" sign is removed the smoothed data will be plotted
    over the unsmoothed data, so that we can she not only the smoothed light curve,
    but also the impact of the smoothing
    '''
    
    fig = plt.figure()
    #plt.plot(data_cleaned[:,0],data_cleaned[:,1],color='orange')
    plt.plot(data_cleaned[2:-2,0],flux_smoothed[2:-2],color='blue')
    plt.xlabel('Time (BJD)')
    plt.ylabel('Flux (electorns/sec)')
    plt.title('Light Curve during Transit')
    plt.legend(['Smoothed data'])
    plt.grid()
    
    if which_txt == 'tess_lc1.dat':
        fig = plt.figure()
        plt.plot(data_cleaned[5720:5810,0],flux_smoothed[5720:5810],color='blue')
        plt.xlabel('Time (BJD)')
        plt.ylabel('Flux (electorns/sec)')
        plt.title('Eclipse Event')
        plt.legend(['Smoothed data'])
        plt.grid()       
    else:
        fig = plt.figure()   
        plt.plot(data_cleaned[6751:7145,0],flux_smoothed[6751:7145],color='blue')
        plt.xlabel('Time (BJD)')
        plt.ylabel('Flux (electorns/sec)')
        plt.title('Eclipse Event')
        plt.legend(['Unsmoothed data', 'Smoothed data'])
        plt.grid()

def eclipse_func(t,f_star, Dt_total, t_cent, Dt_in, d):    
    slope_1 = - d/Dt_in
    slope_2 =  d/Dt_in        
    x_1 = t_cent - 0.5*Dt_total
    x_2 = x_1 + Dt_in
    x_4 = t_cent + 0.5*Dt_total
    x_3 = x_4 - Dt_in

    func = [f_star, lambda t: slope_1 * (t -x_1) + f_star, f_star-d, lambda t: slope_2 * (t -x_4) + f_star, f_star] 
    conditions = [t <= x_1, (x_1 < t) &  (t <= x_2), (x_2 < t) & (t < x_3), (x_3 <= t) & (t < x_4), t >= x_4]
    eclipse_func = np.piecewise(t, conditions, func)
    return eclipse_func
    

#############
### Main ####
#############



### First data sheet

path = os.getcwd() # get the path to working directory

clean_data = data_cleaning(path+'/images/tess_lc1.dat')
smoothed_flux =period_exoplanet(clean_data,path+'/images/tess_lc1.dat')

fit_param_1 , stat = curve_fit(eclipse_func,clean_data[:,0],clean_data[:,1], p0=[178230, 0.1, 0.61+1.333e3, 0.03, 90])
fig = plt.figure()
plt.plot(clean_data[5720:5810,0], clean_data[5720:5810,1])  # 5720:5810 P=1 day
plt.plot(clean_data[5720:5810,0], eclipse_func(clean_data[5720:5810,0],*fit_param_1))
plt.xlabel('Time (BJD)')
plt.ylabel('Flux (electorns/sec)')
plt.title('Fit of the Model on the Data')
plt.legend(['Unsmoothed data','Model'])
plt.grid() 

### Second data sheet    
    
clean_data = data_cleaning(path+'/images/tess_lc2.dat')
smoothed_flux =period_exoplanet(clean_data,path+'/images/tess_lc2.dat')       

fit_param_2 , stat = curve_fit(eclipse_func,clean_data[:,0],clean_data[:,1], p0=[135900, 0.28, 0.3+1.335e3, 0.05, 480])
fig = plt.figure()
plt.plot(clean_data[6751:7145,0], clean_data[6751:7145,1])
plt.plot(clean_data[6751:7145,0], eclipse_func(clean_data[6751:7145,0],*fit_param_2))
plt.xlabel('Time (BJD)')
plt.ylabel('Flux (electorns/sec)')
plt.title('Fit of the Model on the Data')
plt.legend(['Unsmoothed data','Model'])
plt.grid()
     