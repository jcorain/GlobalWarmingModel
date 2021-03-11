# import the different libraries

import numpy as np
import matplotlib.pyplot as plt 

# create a function to define the data frame 

def TemperatureDF(TimeStep: float, 
                WaterDepth: int, 
                L: int, 
                albedo: float, 
                epsilon: float, 
                nSteps: int,
                TempIni: float):
    '''Returns a data frame of the time, temperature and the IR heat from the planet
    

    Parameters
    ----------
    TimeStep : float
        Time step you want to take in years.
    WaterDepth : integer
        The mediunm depth of the Water in m.
    L : integer
        Solar energy constant in W.m-2.
    albedo : float
        albedo reflecting planet energy reflection due to atmosphere layers. Float between 0 and 1
    epsilon : float
        emissivity reflecting the energy emmission of a body regarding the emission of a black body. Float between 0 and 1.
    nSteps : integer
        The number of time steps you are taking.
    TempIni : integer
        The initial temperature

    Returns
    -------
    Numpy array with time, temperature and IR heat form the planet.

    '''
    
    # define stefan boltzman constant
    
    sigma = 5.67E-8 # W.m-2.K-4
    
    # intitiate the different columns of the data frame 
    Time = [0,]
    Temp = [TempIni,]
    HeatOut = [epsilon*sigma*Temp[-1]**4,]
    
    
    # define the palnet heat capacity and the incoming heat 
    
    HeatCapacity = 4200000*WaterDepth
    HeatIn = L*(1-albedo)/4
    HeatContent = HeatCapacity * Temp[0]
    
    # perform the loop for the different time steps
    
    for i in range(nSteps):
        
        Time = np.append(Time, TimeStep*(i+1))
        HeatContent = HeatContent + (HeatIn - HeatOut[-1]) * TimeStep * 3600*24*365
        Temp = np.append(Temp, HeatContent/HeatCapacity)
        HeatOut = np.append(HeatOut, epsilon*sigma*Temp[-1]**4)
          
    # create the data frame at the end and transpose it to have it in a regular way 
    TemperatureArray = np.array([Time,Temp,HeatOut]).transpose()
    return(TemperatureArray)

def plotTemp(TempArray):
    '''Returns a plot of temperature versus time 

    Parameters
    ----------
    TempArray : numpy array
       Numpy array with time, temperature and IR heat form the planet. Can be obtained via the TemperatureDF function

    Returns
    -------
    Plot of temperature versus time.

    '''
    plt.plot(TempArray[:,0], TempArray[:,1])
    plt.xlabel("Time (years)")
    plt.ylabel("Temperature (K)")
    plt.show()
            
        
if __name__ == "__main__":
    # caller for when the script is started
    
    # definne the different parameters
    
    TimeStep = 10 # years
    WaterDepth = 4000 # meters
    L = 1350 # W.m-2
    albedo = 0.3
    epsilon = 1
    nSteps = 200
    TempIni = 0
    
    res = TemperatureDF(TimeStep, 
                      WaterDepth, 
                      L, 
                      albedo, 
                      epsilon,
                      nSteps,
                      TempIni)
    
    # plot the temperature versus time 
    
    plotTemp(res)
    
    # print the latest results 
    
    print(res[-1,1], res[-1,2])
    
     