
from numpy import *  #



def pole(fiA, lamA, fiD, lamD, a, e2):
    b = a * sqrt(1 - e2)                                       #metry
    #print(b)
    fiA = deg2rad(fiA)
    lamA = deg2rad(lamA)
    fiD = deg2rad(fiD)
    lamD = deg2rad(lamD)
    phiA = (sin(fiA) / (1 - e2 * sin(fiA) ** 2)) + (1 / (2 * sqrt(e2))) * log((1 + (sqrt(e2) * sin(fiA))) / (1 - (sqrt(e2) * sin(fiA))))
    #print(phiA)
    phiD = (sin(fiD) / (1 - e2 * sin(fiD) ** 2)) + (1 / (2 * sqrt(e2))) * log((1 + (sqrt(e2) * sin(fiD))) / (1 - (sqrt(e2) * sin(fiD))))
    #print(phiD)

    p = (((b ** 2) * (lamD - lamA)) / 2) * (phiA - phiD)
    return p




def vincent(fiA, lambdaA, fiD, lambdaD, a, e2):
    b = a * sqrt(1 - e2)                                                          # krótsza półoś elipsoidy      #metry
    f = 1 - (b / a)                                                                  # spłaszczenie elipsoidy       #kiczba

    delLambda = (lambdaD - lambdaA)                  # radians(lambdaD - lambdaA)                #różnica długości geodezyjna  #stopni
    Ua = arctan((1 - f) * tan(fiA))                                            # szerokość zredukowana        #stopni
    Ud = arctan((1 - f) * tan(fiD))

    L = delLambda
    while True:
        SIN = sqrt(((cos(Ud) * sin(L)) * (cos(Ud) * sin(L))) +
                         ((cos(Ua) * sin(Ud) - sin(Ua) * cos(Ud) * cos(L))
                          * (cos(Ua) * sin(Ud) - sin(Ua) * cos(Ud) * cos(L))))  # liczba

        COS = sin(Ua) * sin(Ud) + cos(Ua) * cos(Ud) * cos(L)  # liczba
        sigma = arctan(SIN / COS)  # stopnie
        #print(SIN)
        #print(COS)
        #print(sigma)

        SINa = cos(Ua) * cos(Ud) * sin(L) / sin(sigma)  # liczba
        #COS2A = 1 - (SINa * SINa)  # liczba
        COS2A = 1 - (SINa ** 2)  # liczba

        COS2sM = COS - (2 * sin(Ua) * sin(Ud)) / COS2A  # liczba
        C = (f / 16) * COS2A * (4 + f * (4 - 3 * COS2A))  # liczba
        prev_L = L
        L = delLambda + (1 - C) * f * SINa * (sigma + C * SIN * (COS2sM + C * COS * (-1 + 2 * COS2sM)))  # stopnie
        if abs(prev_L - L) < deg2rad(0.000001 / 3600):
            break



    U2 = ((a ** 2 - b ** 2) / b ** 2) * COS2A

    A = 1 + (U2 / 16384) * (4096 + U2 * ((-768) + U2 * (320 - 175 * U2)))

    B = (U2 / 1024) * (256 + U2 * ((-128) + U2 * (74 - 47 * U2)))

    Delsig = B * SIN * (COS2sM + (1 / 4) * B * (COS * ((-1) + 2 * COS2sM ** 2) - (1 / 6) * B * COS2sM * ((-3) + 4 * SIN ** 2) *
                                               ((-3) + 4 * COS2sM ** 2)))

    Sab = b * A * (sigma - Delsig)                                                       #długość linii geodezyjnej pomiędzy punktami A i B

    licznik_ab = cos(Ud) * sin(prev_L)
    mianownik_ab = cos(Ua) * sin(Ud) - sin(Ua) * cos(Ud) * cos(prev_L)

    licznik_ba = cos(Ua) * sin(prev_L)
    #mianownik_ba = sin(Ua) * cos(Ud) + cos(Ua) * sin(Ud) * cos(delLambda)
    mianownik_ba = sin(Ua) * cos(Ud) * cos(prev_L) - sin(Ua) * cos(Ud)

    Aab = arctan(licznik_ab / mianownik_ab)                                       #azymut prosty i odwrotny linii geodezyjnej
    Aba = arctan(licznik_ba / mianownik_ba) + pi


    if mianownik_ab < 0:
        Aab += pi
    elif licznik_ab < 0:
        Aab += 2 * pi
    if Aab > 2 * pi:
        Aab -= 2 * pi
    Aba = Aab + pi


    return Sab, Aab, Aba


Pole = pole(50.25, 20.75, 50.00, 21.25, 6378137, 0.00669437999013)
print(f'Pole powierzchni: {Pole}')

#Vincent = vincent(deg2rad(50.25), deg2rad(20.75), deg2rad(50.00), deg2rad(21.25), 6378137, 0.00669437999013)
#print(Vincent)



def Kivioja(Fi, Lambda, Sab, Azym, a, e2):
    n = int(Sab/1000)
    ds = Sab / n
    for i in range(0, n):
        M = (a * (1 - e2)) / (sqrt((1 - e2 * (sin(Fi) ** 2)) ** 3))
        N = a / sqrt(1 - e2 * (sin(Fi) ** 2))

        DFi = ds * cos(Azym) / M
        FiM = Fi + 0.5 * DFi

        M = (a * (1 - e2)) / sqrt((1 - e2 * (sin(FiM)) ** 2) ** 3)
        N = a / sqrt(1 - e2 * (sin(FiM)) ** 2)

        azM = Azym + ds * sin(Azym) * tan(FiM) / N

        DFi = ds*cos(azM) / M
        Fi += DFi

        DLambda = ds * sin(azM) / (N*cos(FiM))
        Lambda += DLambda
        Daz = ds * sin(azM) * tan(FiM) / N
        Azym += Daz

    return Fi, Lambda, Azym+pi

if __name__ == '__main__':
    # Dane wierzchołków
    fiA = deg2rad(50.25)
    lambdaA = deg2rad(20.75)

    fiB = deg2rad(50.00)
    lambdaB = deg2rad(20.75)

    fiC = deg2rad(50.25)
    lambdaC = deg2rad(21.15)

    fiD = deg2rad(50.00)
    lambdaD = deg2rad(21.25)

    a = 6378137                            # metry             #dłuższa półoś elipsoidy
    e2 = 0.00669437999013                  # bez jednostek

    fiS = (fiA + fiD) / 2
    #print(fiS)
    lambdaS = (lambdaA + lambdaD) / 2
    #print(lambdaS)

    # dane:
    # fi_A, lam_A, fi_D, lam_D
    s_AD, az_AD, az_DA = vincent(fiA, lambdaA, fiD, lambdaD, a, e2)

    fi_S, lam_S, az_SA = Kivioja(fiA, lambdaA, s_AD / 2, az_AD, a, e2)
    SabSR, AabSR, AbaSR = vincent(fiS, lambdaS, fi_S, lam_S, a, e2)

    print(f'odległość: {s_AD} \nAz wprost(Vincenta): {rad2deg(az_AD)} \nAz odwrotny(Vincenta): {rad2deg(az_DA)}')

    print(f'pkt środkowy: {rad2deg(fi_S)} lambda: {rad2deg(lam_S)}')

    print(f'odległość pomiędzy okt środkowymi: {SabSR} \nAz dwoch punktów: {rad2deg(AabSR)}, {rad2deg(AbaSR)}')



