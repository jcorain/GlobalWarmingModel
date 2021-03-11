# import the different libraries

import numpy as np
import matplotlib.pyplot as plt 

# create a funcion that derived the fit parameter fronm albedo and ice latitude regarding temperature 

def FitParam(MeanTempVal):
    '''Returns fit parameter for ice altititde and albedo
    

    Parameters
    ----------
    MeanTempVal : numpy array
        numpy array with three colums : mean temp, ice altitutde and albedo.

    Returns
    -------
    Fit parameter values.

    '''
    
    # get the params for the Ice altitude 
    
    m1, b1 = np.polyfit(MeanTempVal[:,0], MeanTempVal[:,1], 1)
    
    # get the params for the Albedo
    
    m2, b2 = np.polyfit(MeanTempVal[:,0], MeanTempVal[:,2], 1)
    
    return(m1,b1,m2,b2)

# create a function to define the data frame 

def TemperatureAlbedosDF(LRange,
                         LStep: int,
                         albedoIni: float, 
                         epsilon: float,
                         MeanTempVal,
                         nIters: int):
    '''Returns a data frame of the solar Energy, temperature and albedo
    

    Parameters
    ----------
    Lrange : array of two integers
        Array containing max and in Solar energy constant in W.m-2.
    LStep : integer
        The step of solar energy constant you want to take.
    albedoIni : float
        initial albedo reflecting planet energy reflection due to atmosphere layers. Float between 0 and 1
    epsilon : float
        emissivity reflecting the energy emmission of a body regarding the emission of a black body. Float between 0 and 1.
    MeanTempVal : numpy array
        numpy array with three colums : mean temp, ice altitutde and albedo.
    nIters : int
        number of iterations for the inner loop

    Returns
    -------
    Numpy array with time, temperature and IR heat form the planet.

    '''
    
    # define constant
    
    sigma = 5.67E-8 # W.m-2.K-4
    albedoMax = 0.65
    albedoMin = 0.15
    LatMax = 90
    LatMin = 0

    # get the fit values for the albedo versus temperature 
    
    mIce, bIce, mAlbedo, bAlbedo = FitParam(MeanTempVal)
    
    # intitiate the different columns of the data frame 
    
    LArr = [LRange[1],] # Array containing L value
    nIter = [0,] # array containing nIter Value
    albedo = [albedoIni,] # Array containing Albedo
    Temp = [(LArr[0]*(1-albedo[0])/(4*sigma*epsilon))**(1/4),] #Array containing temperature
    IceLat = [mIce*Temp[0]+bIce] #Array containing ice latittude
    HystWay = [-1,] # Array containgn hyterisis way (-1 descending, 1 ascending)

    # create the first loop over L to get the descending values 
    
    # intiate value for L
    
    L = LArr[0]
    
    # do the L decreasing loop 
    
    while L > LRange[0] -1:
        
        if L != LArr[0]:
            # append the latest value to get iter number 0. I use -2 to not take into account the NaN values 
        
            LArr = np.append(LArr, L)
            nIter = np.append(nIter, 0)
            Temp = np.append(Temp, Temp[-2])
            albedo = np.append(albedo, albedo[-2])
            IceLat = np.append(IceLat, IceLat[-2])
            HystWay = np.append(HystWay, -1)
            
        
        # iterate over nIters 
        
        for i in range(nIters):
            
            # append the differnt data to the DF for Larr and nIter
            
            LArr = np.append(LArr, L)
            nIter = np.append(nIter, i+1)
            
            # calulate the temperature with previous value of albedo
            
            Temp = np.append(Temp,
                             (L*(1-albedo[-1])/(4*sigma*epsilon))**(1/4))
            
            # calculate new albedo and Icelat val with last temp value 
            
            albedo = np.append(albedo,
                               min(max(mAlbedo * Temp[-1] + bAlbedo,
                                       albedoMin),
                                   albedoMax))
            IceLat = np.append(IceLat,
                               min(max(mIce*Temp[-1] + bIce,
                                       LatMin),
                                   LatMax))
            
            HystWay = np.append(HystWay, -1)
            
        # add a line with NaN values 
        
        LArr = np.append(LArr, L)
        nIter = np.append(nIter, np.nan)
        Temp = np.append(Temp, np.nan)
        albedo = np.append(albedo, np.nan)
        IceLat = np.append(IceLat, np.nan)
        HystWay = np.append(HystWay, -1)
            
        # decrease L value 
        
        L = L - LStep
        

    # intiate value for L
    
    L = LArr[-2]
    
    # do the L increasing loop 
    
    while L < LRange[1] + 1:
        
        # append the latest value to get iter number 0. I use -2 to not take into account the NaN values 
        
        LArr = np.append(LArr, L)
        nIter = np.append(nIter, 0)
        Temp = np.append(Temp, Temp[-2])
        albedo = np.append(albedo, albedo[-2])
        IceLat = np.append(IceLat, IceLat[-2])
        HystWay = np.append(HystWay, 1)
            
        
        # iterate over nIters 
        
        for i in range(nIters):
            
            # append the differnt data to the DF for Larr and nIter
            
            LArr = np.append(LArr, L)
            nIter = np.append(nIter, i+1)
            
            # calulate the temperature with previous value of albedo
            
            Temp = np.append(Temp,
                             (L*(1-albedo[-1])/(4*sigma*epsilon))**(1/4))
            
            # calculate new albedo and Icelat val with last temp value 
            
            albedo = np.append(albedo,
                               min(max(mAlbedo * Temp[-1] + bAlbedo,
                                       albedoMin),
                                   albedoMax))
            IceLat = np.append(IceLat,
                               min(max(mIce*Temp[-1] + bIce,
                                       LatMin),
                                   LatMax))
            
            HystWay = np.append(HystWay, 1)
            
        # add a line with NaN values 
        
        LArr = np.append(LArr, L)
        nIter = np.append(nIter, np.nan)
        Temp = np.append(Temp, np.nan)
        albedo = np.append(albedo, np.nan)
        IceLat = np.append(IceLat, np.nan)
        HystWay = np.append(HystWay, 1)
            
        # increase L value 
        
        L = L + LStep
    
    return(np.array([LArr,nIter,Temp,albedo,IceLat, HystWay]).transpose())
    

def plotTempAlbedo(TempAlbedoArray,
                    Parameter = "Temperature",
                    plotIter = False,
                    LWay = "asc"
                    ):
    '''Returns a plot of temperature versus time 

    Parameters
    ----------
    TempAlbedoArray : numpy array
        Numpy array obtained via the TemperatureAlbedosDF function
    Parameter : string
        Parameter to plot, either Temperature, Albedo or Ice Latitude
    plotIter : boolean
        Should we plot the differnt iteration or just the hsyterisis form the converged values ?
    LWay : string
        If plotIter, which way to iterate over L, ascending (asc) or descencidng (desc)
    Returns
    -------
    Plots according to options.

    '''
    
    # initiate the lists 
    
    x_list = []
    y_list = []
    
    # restrain the array regarding the LValue
    
    if plotIter == False:
        # get only the latest values form the data 
        DF = TempAlbedoArray[TempAlbedoArray[:,1] == max(TempAlbedoArray[:,1])]
        x_list = DF[:,0]
        xlab = "L (W/m2)"
    elif plotIter == True:
        if LWay == "asc":
            DF = TempAlbedoArray[TempAlbedoArray[:,5] == 1]
        elif LWay == "desc":
            DF = TempAlbedoArray[TempAlbedoArray[:,5] == -1]
        else: 
            print("The LWay is not the expected one.")
            return(None)
        xlab = "Number of iterations"
        x_list = DF[:,1]
    else:
        print("The plotIter parameter is not a boolean.")
        return(None)
        
    
    # define the y_list 
    
    if Parameter == "Temperature":
        y_list = DF[:,2]
        ylab = "Temperature (K)"
    elif Parameter == "Albedo":
        y_list = DF[:,3]
        ylab = "Albedo"
    elif Parameter == "IceLatitude":
        y_list = DF[:,4]
        ylab = "Ice Latitude (° C)"
    else: 
        print("The Parameter does not have the good definition")
        return(None)
    
    # plot the data
    
    plt.plot(x_list, y_list)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.show()
            
        
if __name__ == "__main__":
    # caller for when the script is started
    
    # define the data frame for the Mean temp, Ice Altitiude and Albedo to fit 
    
    MeanTempVal = np.array([[265,255,245,235,225,215],
                           [75,60,45,30,15,0],
                           [0.15,0.25,0.35,0.45,0.55,0.65]]).transpose()
    
    # definne the different parameters
    
    LRange = [1200, 1600] # W.m-2
    LStep = 10 # W.m-2
    albedo = 0.3
    epsilon = 1
    nIters = 100
    
    # L, albedo, nIters = input("").split()
    # L, albedo, nIters = [ float(L), float(albedo), int(nIters) ]
    
    res = TemperatureAlbedosDF(LRange = LRange, 
                               LStep = LStep, 
                               albedoIni = albedo, 
                               epsilon = epsilon,
                               nIters = nIters,
                               MeanTempVal = MeanTempVal)
    
    # plot the temperature versus time 
    
    Parameter = "Temperature"
    plotIter = True
    LWay = "asc"
    
    plotTempAlbedo(TempAlbedoArray=res,
                    Parameter=Parameter,
                    plotIter=plotIter,
                    LWay= LWay)  
    
    # print the latest results 
    
    # restrain the data 
    
    # restDF = res[res[:,0] == L]
     
    # if(nIters == 1):
    #     restDF = restDF[0,:]
    # else:
    #     if albedo == .15:
    #         restDF = restDF[restDF[:,5] == -1]
    #     else:
    #         restDF = restDF[restDF[:,5] == 1]
    #     restDF = restDF[-2,:]
  
    # print(restDF[2], restDF[3])
    
    
  
    
     
