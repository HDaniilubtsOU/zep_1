import math

import numpy as np  #
import pandas as pd  #
import matplotlib.pyplot as plt  #
from matplotlib import cm
from mpl_toolkits import mplot3d

fi = 10

rek = 10.1505555556                                                                           # godziny   (fi)

dek = -0.475833667                                                                           # kąty     (delta)

lam = 20.926389222                                                                           # рассположение обсервотории   (lambda)
h = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]  # czas
year = 2002
mans = 12
day = 20
r = 100

x_list = []
y_list = []
z_list = []


def JD(year, mans, day):
    jd = math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (mans + 1)) + day         # считает от 0.0.0 год
    return jd



def G(Juldate):
    T = (Juldate - 2451545) / 36525
    g = 280.46061837 + 360.98564736629 * (Juldate - 2451545.0) + 0.000387933 * (T ** 2) - (T ** (3/38710000))
    return g


def T(x, GG, lam, rek):                                                           #kąt godzinny
    UT1 = x * 1.002737909350795                                                  # godziny (czas uniwersalny)
    S = UT1 * 15 + lam + GG
    t = S - rek * 15
    return t



    #Rozwiązanie trójkąta paralaktycznego

def ZZ(fi, dek, Tt):
    cosZ = (np.sin(np.radians(fi)) * np.sin(np.deg2rad(dek))) + (np.cos(np.radians(fi)) * np.cos(np.deg2rad(dek)) * np.cos(np.radians(Tt)))       # Rozwiązanie trójkąta paralaktycznego
    Z = np.arccos(cosZ)
    return Z

def AZ(fi, dek, Tt):
    tgA = (-np.cos(np.radians(fi)) * np.sin(np.radians(Tt))) / ((np.cos(np.radians(fi)) * np.sin(np.deg2rad(dek))) -                              #zwraca Azymut w stopniach
                                                                    (np.sin(np.radians(fi)) * np.cos(np.deg2rad(dek)) * np.cos(np.radians(Tt))))
    Az = np.degrees(np.arctan(tgA))

    if (((np.cos(np.radians(fi)) * np.sin(np.deg2rad(dek))) -
             (np.sin(np.radians(fi)) * np.cos(np.deg2rad(dek)) * np.cos(np.radians(Tt)))) < 0):
        Az += 180
    elif ((-np.cos(np.radians(fi)) * np.sin(np.radians(Tt))) < 0):
        Az += 360

    return Az


def Ws(AZz, ZZz):
    #ZZz = np.deg2rad(ZZz)
    #AZz = np.deg2rad(AZz)
    x = r * np.sin((ZZz)) * np.cos(np.deg2rad(AZz))                                                                 # Transformacja współrzędnych
    y = r * np.sin((ZZz)) * np.sin(np.deg2rad(AZz))
    z = r * np.cos((ZZz))
    return x, y, z



for x in h:

    Juldate = JD(2002, 12, 20)
    #print(Juldate)

    GG = G(Juldate)
    #print(GG)

    Tt = T(x, GG, lam, rek)
    #print(Tt)

    ZZz = ZZ(fi, dek, Tt)
    #print(ZZz)

    AZz = AZ(fi, dek, Tt)
    #print(AZz)

    WSs = Ws(AZz, ZZz)
    print(WSs)


    x_list.append(WSs[0])
    y_list.append(WSs[1])
    z_list.append(WSs[2])
#print(x_list)
#print(y_list)
#print(z_list)


fig = plt.figure()
ax = plt.axes(projection="3d")

ax.set_title("Ruch gwiazdY")

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 96 * np.outer(np.cos(u), np.sin(v))
y = 96 * np.outer(np.sin(u), np.sin(v))
z = 96 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b', rstride=3, cstride=5, cmap='seismic', antialiased=True)
ax.scatter3D(x_list, y_list, z_list, c='black', cmap='cividis', s=100)
plt.show()