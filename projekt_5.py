import math
import numpy as np  #
import matplotlib.pyplot as plt

a = 6378137 # metry
e2 = 0.00669437999013 # bez jednostek

ain = 6378245                                 #Krasowskiego
bin = 6356863                                 #Krasowskiego
e2i = (ain ** 2 - bin ** 2) / (ain ** 2)      #Krasowskiego           0.006693427489819221

kappa = 0.8407728e-6
alfa = -1.786877784465417e-06
beta = -2.5612706773016787e-07
gam = 4.089597325631379e-06

xx = -33.4297
yy = 146.5746
zz = 76.2865




def zam1(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z

Z = zam1(50.25, 20.75, 100, 6378137, 0.00669437999013)
print(f'pkt A(x, y, z) GRS80: {Z}')




def BursyWolf1(x, y, z):

    mat_1 = np.array([x, y, z])

    mat_2 = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])

    mat_3 = np.array([xx, yy, zz])

    Ws = mat_1 + mat_2 @ mat_1 + mat_3

    return Ws

BurW = BursyWolf1(3821511.431829001, 1447841.1655117115, 4880693.944015627)
print(f'pkt A(Krasowskiego): {BurW}')


def Hirviona1(x, y, z, ain, e2i):
    r = np.sqrt((x ** 2) + (y ** 2))

    Fi = np.arctan((z / r) * (1 - e2i) ** -1)

    N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)

    h = r / np.cos(Fi) - N


    Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)


    eps = 0.00005
    epsilon = eps

    while abs(np.deg2rad(Fii - Fi) * 3600) >= epsilon:
        Fi = Fii

        N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)
        #print(N)
        h = r / np.cos(Fi) - N
        #print(h)

        Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)

    lam = np.arctan(y / x)

    return np.rad2deg(Fii), np.rad2deg(lam), h

Hir = Hirviona1(3821488.38631705, 1447964.60777067, 4880775.94239304, 6378245, 0.006693427489819221 )
print(f'pkt A(fi, lambda, h) : {Hir}')





def zam2(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z

Z = zam2(50.00, 20.75, 100, 6378137, 0.00669437999013)
print(f'pkt B(x, y, z) GRS80: {Z}')


def BursyWolf2(x, y, z):

    mat_1 = np.array([x, y, z])

    mat_2 = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])

    mat_3 = np.array([xx, yy, zz])

    Ws = mat_1 + mat_2 @ mat_1 + mat_3

    return Ws

BurW = BursyWolf2(3841468.4576248974, 1455402.2062161348, 4862865.642150784)
print(f'pkt B(Krasowskiego): {BurW}')


def Hirviona2(x, y, z, ain, e2i):
    r = np.sqrt((x ** 2) + (y ** 2))

    Fi = np.arctan((z / r) * (1 - e2i) ** -1)

    N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)

    h = r / np.cos(Fi) - N


    Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)


    eps = 0.00005
    epsilon = eps

    while abs(np.deg2rad(Fii - Fi) * 3600) >= epsilon:
        Fi = Fii

        N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)
        #print(N)
        h = r / np.cos(Fi) - N
        #print(h)

        Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)

    lam = np.arctan(y / x)

    return np.rad2deg(Fii), np.rad2deg(lam), h

Hir = Hirviona2(3841445.45524758, 1455525.60507301, 4862947.63393776, 6378245, 0.006693427489819221 )
print(f'pkt B(fi, lambda, h) : {Hir}')





def zam3(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z

Z = zam3(50.25, 21.15, 100, 6378137, 0.00669437999013)
print(f'pkt C(x, y, z) GRS80: {Z}')


def BursyWolf3(x, y, z):

    mat_1 = np.array([x, y, z])

    mat_2 = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])

    mat_3 = np.array([xx, yy, zz])

    Ws = mat_1 + mat_2 @ mat_1 + mat_3

    return Ws

BurW = BursyWolf3(3811310.5482445396, 1474484.8486810413, 4880693.944015627)
print(f'pkt C(Krasowskiego): {BurW}')


def Hirviona3(x, y, z, ain, e2i):
    r = np.sqrt((x ** 2) + (y ** 2))

    Fi = np.arctan((z / r) * (1 - e2i) ** -1)

    N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)

    h = r / np.cos(Fi) - N


    Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)


    eps = 0.00005
    epsilon = eps

    while abs(np.deg2rad(Fii - Fi) * 3600) >= epsilon:
        Fi = Fii

        N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)
        #print(N)
        h = r / np.cos(Fi) - N
        #print(h)

        Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)

    lam = np.arctan(y / x)

    return np.rad2deg(Fii), np.rad2deg(lam), h

Hir = Hirviona3(3811287.6031179, 1474608.35505879, 4880775.99261476, 6378245, 0.006693427489819221 )
print(f'pkt C(fi, lambda, h) : {Hir}')






def zam4(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z

Z = zam4(50.00, 21.25, 100, 6378137, 0.00669437999013)
print(f'pkt D(x, y, z) GRS80: {Z}')


def BursyWolf4(x, y, z):

    mat_1 = np.array([x, y, z])

    mat_2 = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])

    mat_3 = np.array([xx, yy, zz])

    Ws = mat_1 + mat_2 @ mat_1 + mat_3

    return Ws

BurW = BursyWolf4(3828621.5672599915, 1488869.499821071, 4862865.642150784)
print(f'pkt D(Krasowskiego): {BurW}')


def Hirviona4(x, y, z, ain, e2i):
    r = np.sqrt((x ** 2) + (y ** 2))

    Fi = np.arctan((z / r) * (1 - e2i) ** -1)

    N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)

    h = r / np.cos(Fi) - N


    Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)


    eps = 0.00005
    epsilon = eps

    while abs(np.deg2rad(Fii - Fi) * 3600) >= epsilon:
        Fi = Fii

        N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)
        #print(N)
        h = r / np.cos(Fi) - N
        #print(h)

        Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)

    lam = np.arctan(y / x)

    return np.rad2deg(Fii), np.rad2deg(lam), h

Hir = Hirviona4(3828598.69094911, 1488992.97935494, 4862947.69703016, 6378245, 0.006693427489819221 )
print(f'pkt D(fi, lambda, h) : {Hir}')





def zam5(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z

Z = zam5(50.13, 21.00, 100, 6378137, 0.00669437999013)
print(f'pkt S(x, y, z) GRS80: {Z}')


def BursyWolf5(x, y, z):

    mat_1 = np.array([x, y, z])

    mat_2 = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])

    mat_3 = np.array([xx, yy, zz])

    Ws = mat_1 + mat_2 @ mat_1 + mat_3

    return Ws

BurW = BursyWolf5(3824730.291204837, 1468176.4025040697, 4872147.884370603)
print(f'pkt S(Krasowskiego): {BurW}')


def Hirviona5(x, y, z, ain, e2i):
    r = np.sqrt((x ** 2) + (y ** 2))

    Fi = np.arctan((z / r) * (1 - e2i) ** -1)

    N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)

    h = r / np.cos(Fi) - N


    Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)


    eps = 0.00005
    epsilon = eps

    while abs(np.deg2rad(Fii - Fi) * 3600) >= epsilon:
        Fi = Fii

        N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)
        #print(N)
        h = r / np.cos(Fi) - N
        #print(h)

        Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)

    lam = np.arctan(y / x)

    return np.rad2deg(Fii), np.rad2deg(lam), h

Hir = Hirviona5(3824707.32937327, 1468299.86396727, 4872229.91107486, 6378245, 0.006693427489819221 )
print(f'pkt S(fi, lambda, h) : {Hir}')




def zam6(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z

Z = zam6(50.25, 20.75, 100, 6378137, 0.00669437999013)
print(f'pkt K(x, y, z) GRS80: {Z}')


def BursyWolf6(x, y, z):

    mat_1 = np.array([x, y, z])

    mat_2 = np.array(
        [[kappa, gam, -beta],
         [-gam, kappa, alfa],
         [beta, -alfa, kappa]
         ])

    mat_3 = np.array([xx, yy, zz])

    Ws = mat_1 + mat_2 @ mat_1 + mat_3

    return Ws

BurW = BursyWolf6(3821511.431829001, 1447841.1655117115, 4880693.944015627)
print(f'pkt K(Krasowskiego): {BurW}')


def Hirviona6(x, y, z, ain, e2i):
    r = np.sqrt((x ** 2) + (y ** 2))

    Fi = np.arctan((z / r) * (1 - e2i) ** -1)

    N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)

    h = r / np.cos(Fi) - N


    Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)


    eps = 0.00005
    epsilon = eps

    while abs(np.deg2rad(Fii - Fi) * 3600) >= epsilon:
        Fi = Fii

        N = ain / np.sqrt(1 - e2i * np.sin(Fi) ** 2)
        #print(N)
        h = r / np.cos(Fi) - N
        #print(h)

        Fii = np.arctan((z / r) * (1 - e2i * N / (N + h)) ** -1)

    lam = np.arctan(y / x)

    return np.rad2deg(Fii), np.rad2deg(lam), h

Hir = Hirviona6(3821488.38631705, 1447964.60777067, 4880775.94239304, 6378245, 0.006693427489819221 )
print(f'pkt K(fi, lambda, h) : {Hir}')