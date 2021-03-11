# import libraries

import numpy as np
import matplotlib.pyplot as plt
import math
plt.style.use('classic')

def NearFutureDF(yearBeg:int,
                 yearEnd:int,
                 yearStep:int,
                 CO2Eq:int,
                 CO2Beg:int,
                 A:float,
                 drawDown:float,
                 climateSensitivityK:float,
                 climateSensitivity2X:float,
                 RFForcingToday:float,
                 yearNow:int,
                 tResponseTime:int):
    
    # initialize the different columns of the DF 
    
    years = [yearBeg,]
    CO2_BAS = [CO2Beg,]
    RFCO2_BAS = [0,]
    CO2Rate_BAS = [0,]
    RFMasked_BAS = [0,]
    RFTot_BAS = [0,]
    TEq_BAS = [0,]
    TEvol_BAS = [0,]
    CO2_NH = [CO2Beg,]
    RFCO2_NH = [0,]
    RFMasked_NH = [0,]
    RFTot_NH = [0,]
    TEq_NH = [0.]
    TEvol_NH = [0,]
    
    # create a loop over years and calulate the CO2 raise in the business as usual scenario
    
    for i in range(int((yearEnd-yearBeg)/yearStep)):
        years = np.append(years, years[-1] + yearStep)
        CO2_BAS = np.append(CO2_BAS, CO2Eq + (CO2_BAS[-1] - CO2Eq)*(1+A*yearStep))
        RFCO2_BAS = np.append(RFCO2_BAS, climateSensitivity2X*math.log(CO2_BAS[i]/CO2Eq)/math.log(2))
        CO2Rate_BAS = np.append(CO2Rate_BAS, (CO2_BAS[-1]-CO2_BAS[-2])/yearStep)

    
    # calculate B 
    
    idxToday = np.where(years== yearNow)[0]
    B = RFForcingToday/CO2Rate_BAS[idxToday]
    
    # continue the time loop using the B value 
    
    for i in range(int((yearEnd-yearBeg)/yearStep)):
        RFMasked_BAS = np.append(RFMasked_BAS, max(B*CO2Rate_BAS[i+1],RFForcingToday))
        RFTot_BAS = np.append(RFTot_BAS, RFCO2_BAS[i+1] + RFMasked_BAS[i+1])
        TEq_BAS = np.append(TEq_BAS, climateSensitivityK/climateSensitivity2X*RFTot_BAS[i+1])
        TEvol_BAS = np.append(TEvol_BAS,TEvol_BAS[-1]+(TEq_BAS[-1]-TEvol_BAS[-1])*yearStep/tResponseTime)
        
        # create the loop for the non human scenario 
        
        if i < idxToday:
            CO2_NH = np.append(CO2_NH, CO2_BAS[i+1])
            RFCO2_NH = np.append(RFCO2_NH, RFCO2_BAS[i+1])
            RFMasked_NH = np.append(RFMasked_NH, RFMasked_BAS[i+1])
            RFTot_NH = np.append(RFTot_NH, RFTot_BAS[i+1])
        else:
            CO2_NH = np.append(CO2_NH, CO2_NH[-1] + (CO2Eq - CO2_NH[-1])*drawDown*yearStep)
            RFCO2_NH = np.append(RFCO2_NH, climateSensitivity2X*math.log(CO2_NH[i+1]/CO2Eq)/math.log(2))
            RFMasked_NH = np.append(RFMasked_NH, 0)
            RFTot_NH = np.append(RFTot_NH, RFCO2_NH[i+1])
        
        # calulate the different temperature 
        
        TEq_NH = np.append(TEq_NH, climateSensitivityK/climateSensitivity2X*RFTot_NH[i+1])
        TEvol_NH = np.append(TEvol_NH,TEvol_NH[-1]+(TEq_NH[-1]-TEvol_NH[-1])*yearStep/tResponseTime)
        
    TDiff = TEvol_BAS - TEvol_NH

    out = np.vstack((years, CO2_BAS, RFCO2_BAS, CO2Rate_BAS, RFMasked_BAS, RFTot_BAS, TEq_BAS, TEvol_BAS,
                     CO2_NH, RFCO2_NH, RFMasked_NH, RFTot_NH, TEq_NH, TEvol_NH, TDiff)).transpose()
    
    
        
    return(out)

def NearFuturePlot(DF,
                   yParam:str):
    
    # define the x valuesas years
    
    x = DF[:,0]
    
    # define the y values regarding yparam
    
    if yParam == "CO2":
        y1 = DF[:,1]
        y2 = DF[:,8]
        ylab = "C02 emission (ppm)"
        legend = ["Business as usual","Human emission stopped at " + str(yearNow)]
    elif yParam == "RF":
        y1 = DF[:,5]
        y2 = DF[:,4]
        y3 = DF[:,11]
        y4 = DF[:,10]
        ylab = "Radiative forcing (W.m-2)"
        legend = ["Business as usual: Total","Business as usual: Masked",
                  "Human emission stopped at " + str(yearNow) + ": Total","Human emission stopped at " + str(yearNow) + ": Masked"]
    elif yParam == "TemperatureRaise":
        y1 = DF[:,7]
        y2 = DF[:,13]
        ylab = "Temperature Raise (K)"
        legend = ["Business as usual","Human emission stopped at " + str(yearNow)]
        
        
    # plot the data 
    
    if yParam == "RF":
        plt.plot(x, y1, x, y2, x, y3, x, y4)
    else:
        plt.plot(x, y1, x, y2)
    plt.xlabel("Years")
    plt.ylabel(ylab)
    plt.legend(legend)
    
if __name__ == "__main__":
    
    # define the constant 
    
    yearBeg = 1900 # years
    yearEnd = 2100 #years
    yearStep = 1 #years
    CO2Eq = 280 #ppm
    CO2Beg = 290 #ppm
    A = 0.0225 # year-1
    drawDown = 0.01 #year-1
    climateSensitivityK = 6.15 # K
    climateSensitivity2X = 4 # Wm-2
    RFForcingToday = -1.5 # Wm-2
    yearNow = 2015 # year
    tResponseTime = 20 # years 
    yParam = "TemperatureRaise"
  
    res = NearFutureDF(yearBeg=yearBeg, 
                       yearEnd=yearEnd, 
                       yearStep = yearStep, 
                       CO2Eq = CO2Eq, 
                       CO2Beg = CO2Beg, 
                       A = A, 
                       drawDown=drawDown, 
                       climateSensitivityK=climateSensitivityK,
                       climateSensitivity2X=climateSensitivity2X,
                       RFForcingToday=RFForcingToday,
                       yearNow=yearNow,
                       tResponseTime=tResponseTime)
    
    NearFuturePlot(DF=res, 
                   yParam=yParam)
    
    
    
    yearIdx = np.where(res[:,0] == yearNow)[0]
    
    print(res[yearIdx, 7])