from numpy import *  #



a = 6378137 # metry
e2 = 0.00669437999013 # bez jednostek

# Dane wierzchołków
fiA = deg2rad(50.25)
#print(fiA)
lambdaA = deg2rad(20.75)

fiB = deg2rad(50.00)
lambdaB = deg2rad(20.75)

fiC = deg2rad(50.25)
lambdaC = deg2rad(21.15)

fiD = deg2rad(50.00)
lambdaD = deg2rad(21.25)

fiSr = deg2rad(50.12525817437701)
lambdaSr = deg2rad(21.000636239332643)

fiS = (fiA + fiD) / 2
lamS = (lambdaA + lambdaD) / 2

f = 1 / 298.257222101
b = ( a - ( f * a ) )

ep2 = ((a ** 2) - (b ** 2)) / (b ** 2)
#print(ep2)


def pole(fiA, lamA, fiD, lamD, a, e2):
    b = a * sqrt(1 - e2)                                       #metry
    #print(b)
    phiA = (sin(fiA) / (1 - e2 * sin(fiA) ** 2)) + (1 / (2 * sqrt(e2))) * log((1 + (sqrt(e2) * sin(fiA))) / (1 - (sqrt(e2) * sin(fiA))))
    #print(phiA)
    phiD = (sin(fiD) / (1 - e2 * sin(fiD) ** 2)) + (1 / (2 * sqrt(e2))) * log((1 + (sqrt(e2) * sin(fiD))) / (1 - (sqrt(e2) * sin(fiD))))
    #print(phiD)

    p = (((b ** 2) * (lamD - lamA)) / 2) * (phiA - phiD)
    return p

Pole = pole(50.25, 20.75, 50.00, 21.25, 6378137, 0.00669437999013)
#print(Pole)




def liczM(fi):
    return a * (1 - e2) / (sqrt(1 - e2 * (sin(fi)) ** 2)) ** 3


def liczN(fi):
    return a / sqrt(1 - e2 * (sin(fi)) ** 2)


def GKruger(fi, la, L0):
    L0 *= pi / 180

    A0 = 1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256
    A2 = 3 * (e2 + e2 ** 2 / 4 + 15 * e2 ** 3 / 128) / 8
    A4 = 15 * (e2 ** 2 + 3 * e2 ** 3 / 4) / 256
    A6 = 35 * e2 ** 3 / 3072
    o = a * (A0 * fi - A2 * sin(2 * fi) + A4 * sin(4 * fi) - A6 * sin(6 * fi))

    t = tan(fi)
    n2 = ep2 * cos(fi) ** 2
    l = la - L0
    Xgk = o + 0.5 * l ** 2 * liczN(fi) * sin(fi) * cos(fi) * (
                1 + l ** 2 / 12 * cos(fi) ** 2 * (5 - t ** 2 + 9 * n2 + 4 * n2 ** 2) + l ** 4 / 360 * cos(fi) ** 4 * (
                    61 - 58 * t ** 2 + t ** 4 + 270 * n2 - 330 * n2 * t ** 2))
    Ygk = l * liczN(fi) * cos(fi) * (1 + l ** 2 / 6 * cos(fi) ** 2 * (1 - t ** 2 + n2) + l ** 4 / 120 * cos(fi) ** 4 * (
                5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2))
    #print(n2)

    return Xgk, Ygk


def FL92(fi, la, acc=3):
    Xgk, Ygk = GKruger(fi, la, 19)
    x92 = 0.9993 * Xgk - 5300000
    y92 = 0.9993 * Ygk + 500000
    x92 = round(x92, acc)
    y92 = round(y92, acc)

    return x92, y92


def strefa(la):
    la *= 180 / pi
    la += 0.5
    la = int(la)
    return la / 3


def FL2k(fi, la, acc=3):
    nrS = strefa(la)
    Xgk, Ygk = GKruger(fi, la, nrS * 3)
    x2k = 0.999923 * Xgk
    y2k = 0.999923 * Ygk + 500000 + nrS * 1000000
    x2k = round(x2k, acc)
    y2k = round(y2k, acc)

    return x2k, y2k


def GKFL(x, y, L0):
    A0 = 1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256
    A2 = 3 * (e2 + e2 ** 2 / 4 + 15 * e2 ** 3 / 128) / 8
    A4 = 15 * (e2 ** 2 + 3 * e2 ** 3 / 4) / 256
    A6 = 35 * e2 ** 3 / 3072

    mianow = a * A0
    fi1 = x / mianow
    sigma = a * (A0 * fi1 - A2 * sin(2 * fi1) + A4 * sin(4 * fi1) - A6 * sin(6 * fi1))
    o = sigma

    eps = deg2rad(0.000001 / 3600)
    while True:
        old = fi1
        fi1 += (x - o) / mianow
        o = sigma

        if abs(old - fi1) < eps:
            break

    N = liczN(fi1)
    M = liczM(fi1)
    t = tan(fi1)
    n2 = ep2 * cos(fi1) ** 2

    fi = fi1 - y ** 2 * t / (2 * M * N) * (
                1 - y ** 2 / (12 * N ** 2) * (5 + 3 * t ** 2 + n2 - 9 * n2 * t ** 2 - 4 * n2 ** 2) + y ** 4 / (
                    360 * N ** 4) * (61 + 90 * t ** 2 + 45 * t ** 4))
    la = L0 * pi / 180 + y / (N * cos(fi1)) * (
                1 - y ** 2 / (6 * N ** 2) * (1 + 2 * t ** 2 + n2) + y ** 4 / (120 * N ** 4) * (
                    5 + 28 * t ** 2 + 24 * t ** 4 + 6 * n2 + 8 * n2 * t ** 2))


    return fi, la


def u92FL(x, y):
    Xgk = (x + 5300000) / 0.9993
    Ygk = (y - 500000) / 0.9993
    return GKFL(Xgk, Ygk, 19)


def u2kFL(x, y):
    Xgk = x / 0.999923
    L0 = int((y - 500000) / 1000000) * 3
    Ygk = (y - 500000) % 1000000 / 0.999923
    return GKFL(Xgk, Ygk, L0)



def polexy(x1, y1, x2, y2):
    return abs((x2 - x1) * (y2 - y1))



if __name__ == '__main__':
    xgkA, ygkA = GKruger(fiA, lambdaA, 19)
    xA, yA = FL92(fiA, lambdaA)                                         # układ 92
    x2kA, y2kA = FL2k(fiA, lambdaA)
    #print(ygkA)

    xgkB, ygkB = GKruger(fiB, lambdaB, 19)
    xB, yB = FL92(fiB, lambdaB)  # układ 92
    x2kB, y2kB = FL2k(fiB, lambdaB)

    xgkC, ygkC = GKruger(fiC, lambdaC, 19)
    xC, yC = FL92(fiC, lambdaC)  # układ 92
    x2kC, y2kC = FL2k(fiC, lambdaC)

    xgkD, ygkD = GKruger(fiD, lambdaD, 19)
    xD, yD = FL92(fiD, lambdaD)  # układ 92
    x2kD, y2kD = FL2k(fiD, lambdaD)

    xgkS, ygkS = GKruger(fiS, lambdaSr, 19)
    xS, yS = FL92(fiS, lambdaSr)  # układ 92
    x2kS, y2kS = FL2k(fiS, lambdaSr)

    xgkK, ygkK = GKruger(fiS, lamS, 19)
    xK, yK = FL92(fiS, lamS)  # układ 92
    x2kK, y2kK = FL2k(fiS, lamS)

    R = sqrt(liczM(fiA) * liczN(fiA))
    #print(R)
    mA = 1 + ygkA ** 2 / (2 * R ** 2) + ygkA ** 4 / (24 * R ** 4)
    KA = (mA - 1) * 1000
    m92A = mA * 0.9993
    K92A = (m92A - 1) * 1000
    m2kA = mA * 0.999923
    K2kA = (m2kA - 1) * 1000

    mA2 = mA ** 2
    kA2 = (mA2 - 1) * 10000
    m92A2 = mA2 * 0.9993
    K92A2 = (m92A2 - 1) * 1000
    m2kA2 = mA2 * 0.999923
    K2kA2 = (m2kA2 - 1) * 1000



    R = sqrt(liczM(fiB) * liczN(fiB))
    mB = 1 + ygkB ** 2 / (2 * R ** 2) + ygkB ** 4 / (24 * R ** 4)
    KB = (mB - 1) * 1000
    m92B = mB * 0.9993
    K92B = (m92B - 1) * 1000
    m2kB = mB * 0.999923
    K2kB = (m2kB - 1) * 1000

    mB2 = mB ** 2
    kB2 = (mB2 - 1) * 10000
    m92B2 = mB2 * 0.9993
    K92B2 = (m92B2 - 1) * 1000
    m2kB2 = mB2 * 0.999923
    K2kB2 = (m2kB2 - 1) * 1000
    
    

    R = sqrt(liczM(fiC) * liczN(fiC))
    mC = 1 + ygkC ** 2 / (2 * R ** 2) + ygkC ** 4 / (24 * R ** 4)
    KC = (mC - 1) * 1000
    m92C = mC * 0.9993
    K92C = (m92C - 1) * 1000
    m2kC = mC * 0.999923
    K2kC = (m2kC - 1) * 1000

    mC2 = mC ** 2
    kC2 = (mC2 - 1) * 10000
    m92C2 = mC2 * 0.9993
    K92C2 = (m92C2 - 1) * 1000
    m2kC2 = mC2 * 0.999923
    K2kC2 = (m2kC2 - 1) * 1000
    
    

    R = sqrt(liczM(fiD) * liczN(fiD))
    mD = 1 + ygkD ** 2 / (2 * R ** 2) + ygkD ** 4 / (24 * R ** 4)
    KD = (mD - 1) * 1000
    m92D = mD * 0.9993
    K92D = (m92D - 1) * 1000
    m2kD = mD * 0.999923
    K2kD = (m2kD - 1) * 1000

    mD2 = mD ** 2
    kD2 = (mD2 - 1) * 10000
    m92D2 = mD2 * 0.9993
    K92D2 = (m92D2 - 1) * 1000
    m2kD2 = mD2 * 0.999923
    K2kD2 = (m2kD2 - 1) * 1000
    
    

    R = sqrt(liczM(fiSr) * liczN(fiSr))
    mSr = 1 + ygkS ** 2 / (2 * R ** 2) + ygkS ** 4 / (24 * R ** 4)
    KSr = (mSr - 1) * 1000
    m92Sr = mSr * 0.9993
    K92Sr = (m92Sr - 1) * 1000
    m2kSr = mSr * 0.999923
    K2kSr = (m2kSr - 1) * 1000

    mSr2 = mSr ** 2
    kSr2 = (mSr2 - 1) * 10000
    m92Sr2 = mSr2 * 0.9993
    K92Sr2 = (m92Sr2 - 1) * 1000
    m2kSr2 = mSr2 * 0.999923
    K2kSr2 = (m2kSr2 - 1) * 1000
    
    

    R = sqrt(liczM(fiS) * liczN(fiS))
    mS = 1 + ygkS ** 2 / (2 * R ** 2) + ygkS ** 4 / (24 * R ** 4)
    KS = (mS - 1) * 1000
    m92S = mS * 0.9993
    K92S = (m92S - 1) * 1000
    m2kS = mS * 0.999923
    K2kS = (m2kS - 1) * 1000

    mS2 = mS ** 2
    kS2 = (mS2 - 1) * 10000
    m92S2 = mS2 * 0.9993
    K92S2 = (m92S2 - 1) * 1000
    m2kS2 = mS2 * 0.999923
    K2kS2 = (m2kS2 - 1) * 1000


    print(f'ws. pkA GK: {xgkA}, {ygkA}')
    print(f'ws. pkB GK: {xgkB}, {ygkB}')
    print(f'ws. pkC GK: {xgkC}, {ygkC}')
    print(f'ws. pkD GK: {xgkD}, {ygkD}')
    print(f'ws. pkS GK: {xgkS}, {ygkS}')
    print(f'ws. pkK GK: {xgkK}, {ygkK}')

    print(f'ws. pkA G1992: {xA}, {yA}')
    print(f'ws. pkB G1992: {xB}, {yB}')
    print(f'ws. pkC G1992: {xC}, {yC}')
    print(f'ws. pkD G1992: {xD}, {yD}')
    print(f'ws. pkS G1992: {xS}, {yS}')
    print(f'ws. pkK G1992: {xK}, {yK}')

    print(f'ws. pkA 2000: {x2kA}, {y2kA}')
    print(f'ws. pkB 2000: {x2kB}, {y2kB}')
    print(f'ws. pkC 2000: {x2kC}, {y2kC}')
    print(f'ws. pkD 2000: {x2kD}, {y2kD}')
    print(f'ws. pkS 2000: {x2kS}, {y2kS}')
    print(f'ws. pkK 2000: {x2kK}, {y2kK}')


    poleGK = polexy(xgkA, ygkA, xgkD, ygkD)
    pole1992 = polexy(xA, yA, xD, yD)
    pole2000 = polexy(x2kA, y2kA, x2kD, y2kD)

    print(f'pole GK: {poleGK}')
    print(f'pole 1992: {pole1992}')
    print(f'pole 2000: {pole2000}')

    #print(mA)#, KA, m92A, K92A, m2kA, K2kA)
    print(f'mA: {mA}')
    print(f'mB: {mB}')
    print(f'mC: {mC}')
    print(f'mD: {mD}')
    print(f'mS: {mSr}')
    print(f'mK: {mS}')

    print(f'KA: {KA}')
    print(f'KB: {KB}')
    print(f'KC: {KC}')
    print(f'KD: {KD}')
    print(f'KSr: {KSr}')
    print(f'Ks: {KS}')

    print(f'm92A: {m92A}')
    print(f'm92B: {m92B}')
    print(f'm92C: {m92C}')
    print(f'm92D: {m92D}')
    print(f'm92Sr: {m92Sr}')
    print(f'm92s: {m92S}')

    print(f'K92A: {K92A}')
    print(f'K92B: {K92B}')
    print(f'K92C: {K92C}')
    print(f'K92D: {K92D}')
    print(f'K92Sr: {K92Sr}')
    print(f'K92s: {K92S}')

    print(f'm2kA: {m2kA}')
    print(f'm2kB: {m2kB}')
    print(f'm2kC: {m2kC}')
    print(f'm2kD: {m2kD}')
    print(f'm2kSr: {m2kSr}')
    print(f'm2ks: {m2kS}')

    print(f'K2kA: {K2kA}')
    print(f'K2kB: {K2kB}')
    print(f'K2kC: {K2kC}')
    print(f'K2kD: {K2kD}')
    print(f'K2kSr: {K2kSr}')
    print(f'K2ks: {K2kS}')

    print(f'mA^2: {mA2}')
    print(f'mB^2: {mB2}')
    print(f'mC^2: {mC2}')
    print(f'mD^2: {mD2}')
    print(f'mS^2: {mSr2}')
    print(f'mK^2: {mS2}')

    print(f'KA^2: {kA2}')
    print(f'KB^2: {kB2}')
    print(f'KC^2: {kC2}')
    print(f'KD^2: {kD2}')
    print(f'KSr^2: {kSr2}')
    print(f'Ks^2: {kS2}')

    print(f'm92A^2: {m92A2}')
    print(f'm92B^2: {m92B2}')
    print(f'm92C^2: {m92C2}')
    print(f'm92D^2: {m92D2}')
    print(f'm92Sr^2: {m92Sr2}')
    print(f'm92s^2: {m92S2}')

    print(f'K92A^2: {K92A2}')
    print(f'K92B^2: {K92B2}')
    print(f'K92C^2: {K92C2}')
    print(f'K92D^2: {K92D2}')
    print(f'K92Sr^2: {K92Sr2}')
    print(f'K92s^2: {K92S2}')

    print(f'm2kA^2: {m2kA2}')
    print(f'm2kB^2: {m2kB2}')
    print(f'm2kC^2: {m2kC2}')
    print(f'm2kD^2: {m2kD2}')
    print(f'm2kSr^2: {m2kSr2}')
    print(f'm2ks^2: {m2kS2}')

    print(f'K2kA^2: {K2kA2}')
    print(f'K2kB^2: {K2kB2}')
    print(f'K2kC^2: {K2kC2}')
    print(f'K2kD^2: {K2kD2}')
    print(f'K2kSr^2: {K2kSr2}')
    print(f'K2ks^2: {K2kS2}')



