# -*- coding: utf-8 -*-

import numpy as np
import math
import json
import rasterio
from matplotlib import pyplot

def Llambda(im,ML,AL):
    SpectralRadiance = np.zeros((x,y))
    for i in range(x):
        for j in range(y):
            SpectralRadiance[i,j] = ML * im[i,j] + AL
    return SpectralRadiance

def TB(Llambda, K1, K2):
    BrightnessT = np.zeros((x,y))
    for i in range(x):
        for j in range(y):
            pix = Llambda[i,j]
            BrightnessT[i,j] = K2/(math.log((K1/pix)+1)) -273.15
    return BrightnessT

def calcNDVI(im):
    NIR = im[constant['bandNIR'],:,:]
    RED = im[constant['bandRED'],:,:]
    NDVI = (NIR-RED)/(NIR+RED)
    return NDVI

def calcPv(NDVI):
    VegetationRatio = np.zeros((x,y))
    minNDVI = np.min(NDVI)
    maxNDVI = np.max(NDVI)
    for i in range(x):
        for j in range(y):
            VegetationRatio[i,j] = (NDVI[i,j] - minNDVI)/((maxNDVI - minNDVI)**2)
    return VegetationRatio

def Emi(Pv):
    Emissivite = np.zeros((x,y))
    for i in range(x):
        for j in range(y):
            Emissivite[i,j] = Pv[i,j]*0.004+0.986
    return Emissivite

def T(Tb,emi):
    T = np.zeros((x,y))
    c2 = ((6.3*10**-34)*(3*10**8))/(1.4*10**(-23))
    la = constant['lambda']    
    for i in range(x):
        for j in range(y):
            T[i,j] = Tb[i,j]/(1+(Tb[i,j]/c2)*la*math.log(emi[i,j]))
    return T

def calcHFI(T):
    HFI = np.zeros((x,y))
    minT = np.min(T)
    maxT = np.max(T)
    for i in range(x):
        for j in range(y):
            HFI[i,j] = (T[i,j]-minT)/(maxT-minT)
    return HFI

def HeatLevel(pix):
    if pix>=0 and pix<0.15:
        return 0
    elif pix>=0.15 and pix<0.30:
        return 1
    elif pix>=0.30 and pix<0.45:
        return 2
    elif pix>=0.45 and pix<0.60:
        return 3
    elif pix>=0.60 and pix<0.75:
        return 3
    elif pix>=0.75 and pix<0.90:
        return 4
    elif pix>=0.90 and pix<=1:
        return 5
    else:
        print(pix)

def calcUHI(HFI):
    UHI = np.zeros((x,y))
    for i in range(x):
        for j in range(y):
            UHI[i,j] = HeatLevel(HFI[i,j])
    return UHI

# read the data
with open("constants.json", "r") as constants:
    constant = json.load(constants)

# read the image    
image = rasterio.open("image.tif").read()

# read the thermal band
imageThermal = image[constant['bandThermal'],:,:]

[x,y] = imageThermal.shape


SpectralRadiance = Llambda(imageThermal,constant['ML'],constant['AL'])

BrighnessTemperature = TB(SpectralRadiance, constant['K1'],constant['K2'])

NDVI = calcNDVI(image)

Pv = calcPv(NDVI)

Emissivite = Emi(Pv)

LST = T(BrighnessTemperature,Emissivite)

HFI = calcHFI(LST)

UHI = calcUHI(HFI)

pyplot.imshow(UHI)







































