# import the different libraries

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as ani

def IceSheetDF(nX:int,
               domainWidth:int,
               timeStep:int,
               nYears:int,
               flowParam:float,
               snowFall:float):
    '''Returns DF for the Ice Sheet flow model
    

    Parameters
    ----------
    nX : int
        Number of grid points.
    domainWidth : int
        Width of the Ice Sheet in meters.
    timeStep : int
        Years for the time step.
    nYears : int
        Total number of years to span over.
    flowParam : float
        The flow rate of the ice in meter per year.
    snowFall : float
        The fall rate of the snow in meter per year.

    Returns
    -------
    Data frame containing the elevation values of the ice sheet over years.

    '''

    # define the grid steps  
    
    dX = domainWidth/nX

    # initialize the data
    
    years = [0,]
    flow = [0 for i in range(nX+2)]
    elevation = [0 for i in range(nX+2)]
    elevationArr = elevation
    
    # loop over the years 
    
    while years[-1] < nYears:
        years = np.append(years, years[-1] + timeStep)
        
        # loop over the elevations 
        
        for ix in range(nX+2):
            
            # set the values to 0 for the second ghost val
            
            if ix == nX+1:
                elevation[ix] = 0
                flow[ix] = 0
            else:
                
                # iterate over the flow
                
                flow[ix] = (elevation[ix] - elevation[ix+1])/dX * flowParam*(elevation[ix]+elevation[ix+1])/2/dX
            
                # iterate the elevation value setting the first and last value to 0
            
                if ix == 0:
                    elevation[ix] = 0
                else:
                    elevation[ix] = elevation[ix] + (snowFall + flow[ix-1] - flow[ix])*timeStep
                
        # append elevation to elevation array
        
        elevationArr = np.vstack((elevationArr, elevation))
   
   # bundle years and elevationArr     
   
    out = np.hstack((years[:,np.newaxis],elevationArr))
    
   
    return(out)

def IceSheetPlot(IceSheetDF):
    '''Returns animated plot of the ice sheet development over years
    

    Parameters
    ----------
    IceSheetDF : numpy array
        Data frame containing the different elevations values over the years.

    Returns
    -------
    Plot of elevation over years.

    '''
    
    # initialize the plots 
    
    # automatize plotLimit
    
    # delete NaN 
    
    maxVal = np.amax(IceSheetDF[~np.isnan(IceSheetDF).any(axis = 1)][:,1:])
    
    plotLimit = round(maxVal,-len(str(round(maxVal))) +1) + 10**(len(str(round(maxVal))) - 1)
        
    
    # define the plot 
    
    fig,ax =  plt.subplots()
    

    # define the animation 

    def animate(i):
        ax.clear()
        ax.set_ylim([0,plotLimit])
        ax.set_xlabel("nX")
        ax.set_ylabel("Elevation (m)")
        plt.plot(IceSheetDF[i,1:])
        plt.annotate("Year = " + str(IceSheetDF[i,0]), (0,plotLimit-500))

    animator = ani.FuncAnimation(fig, animate, frames = (IceSheetDF.shape[0]), repeat = False, interval = 1)

    plt.show()
    
   
if __name__ == "__main__":
    
    # define the differnet parameters
    
    nX = 10                # number of grid points
    domainWidth = 2e6      # meters
    timeStep = 50        # years
    nYears = 50000        # years
    flowParam = 1e4        # m horizontal / yr
    snowFall = 1         # m / y
  

    # get the DF
    
    res = IceSheetDF(nX = nX,
                     domainWidth=domainWidth,
                     timeStep=timeStep,
                     nYears=nYears,
                     flowParam=flowParam,
                     snowFall=snowFall)

    # print the results
    
    IceSheetPlot(res) 
    
    # get the last resulst 
    
    print(res[-1,1 + round(nX/2)])