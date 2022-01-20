import math

import numpy as np  #
import pandas as pd  #
import matplotlib.pyplot as plt  #
from mpl_toolkits import mplot3d

lot = np.loadtxt("geo.w.txt") # Boryspil-->London
dane_lotu = pd.DataFrame(lot)
plt.scatter(dane_lotu[1], dane_lotu[0])
plt.show()
#print(lot)


Bu_fi, Bu_lambda, Bu_h = 50.3578,30.9063, 137        # Lotnisko BORYSPIL, UKRAINE
Luk_fi, Luk_lambda, Luk_h = 52.2191, 16.4533, 11582  # Lotnisko LONDON, UNITED KINGDOM

a = 6378137 # metry / параметр для grs80
e2 = 0.00669438002290 # bez jednostek / нужна для xyz / mimośrud


def zam(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5

    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))

    return x, y, z


def geo2neu(f, l, h, ff, ll, hh):
    pkt_1 = zam(f, l, h, a, e2)   # konwertujemy/конвертация
    pkt_2 = zam(ff, ll, hh, a, e2)  # konwertujemy/конвертанця

    f = np.deg2rad(f)   # угол на рад
    l = np.deg2rad(l)
    R = np.array(
                [[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                 [-np.sin(f) * np.sin(l), np.cos(l), np.cos(f) * np.sin(l)],
                 [np.cos(f), 0, np.sin(f)]
                 ])

    R = R.transpose()   # транспонирование матрицы обращения
    x = np.array([[pkt_2[0] - pkt_1[0]],  # p.array делает таблицу/матрицу
                  [pkt_2[1] - pkt_1[1]],
                  [pkt_2[2] - pkt_1[2]]])
    neu = R @ x   # @ = оператор умножения матриц
    return neu    # получаем матрицу


n = []
e = []
u = []

# na neu
for i in lot:
    temp = geo2neu(Bu_fi, Bu_lambda, Bu_h, i[0], i[1], i[2])
    n.append(temp[0][0])
    e.append(temp[1][0])
    u.append(temp[2][0])

tan_a = []
s = []
cos_z = []

plt.title('Graf neu')
plt.plot(range(len(u)), u)

for ij in range(len(n)):
    tan_a.append(math.atan(e[ij] / n[ij])) # значение тангенса/tangens   ???
    #a.append(math.atan)
    Sji = (n[ij]**2 + e[ij]**2 + u[ij]**2)**0.5     # значение S
    s.append(Sji)
    cos_z.append(u[ij] / Sji)

    print(tan_a[ij], s[ij], cos_z[ij])              # кривая аз


def azymut(n, e, func):
    Aij = func(e / n)
    if n > 0 and e > 0:
        Aij = np.rad2deg(Aij)
    elif n < 0 and e > 0:
        Aij = np.rad2deg(Aij + np.pi)
    elif n < 0 and e < 0:
        Aij = np.rad2deg(Aij + np.pi)
    else:
            Aij = np.rad2deg(Aij + 2 * np.pi)
    if Aij > 360:
        Aij -= 360
    elif Aij < 0:
        Aij += 360
    return Aij
print(azymut)



fig = plt.figure()
ax = plt.axes(projection="3d")
ax.scatter3D(n, e, u, c=u, cmap='inferno', edgecolor='none', linewidth=1.5)
ax.set_title('FlightAware (neu)')
ax.legend()


for i in range(0, len(u)):
  if u[i] < 0:
       print("zniknięcie", n[i], e[i], u[i]) # самолёт пропадает в этом месте
       break

plt.show()

