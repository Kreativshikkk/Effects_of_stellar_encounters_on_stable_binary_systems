#!/usr/bin/env python3

from sys import argv
import os.path
#import matplotlib.pyplot as plt
import numpy as np
import random
import shutil
#import ipdb

data = argv
m = np.array([1, 1, 0])
r = np.array([float(data[1]) / 2, float(data[1]) / 2, 0])
m[2] = float(data[2])
x_new = float(data[3])
y_new = float(data[4])
v = float(data[5])
steps = int(data[6])
dt = float(data[7])
angle = random.uniform(0, 2 * np.pi)
G = 2.959669916e-4

#массивы нолей для всех величин (z так и остается нулевым)
x = np.empty(len(m))
y = np.empty(len(m))
z = np.empty(len(m))
vx = np.empty(len(m))
vy = np.empty(len(m))
vz = np.empty(len(m))

filename = data[8][:-1]
print(filename)

#функция сохранения данных
def save(filename, x, y, z, vx, vy, vz, m, t, K, P):
    np.savez(filename, x = x, y = y, z = z, vx = vx, vy = vy, vz = vz, m = m, t = t, K = K, P = P)
    return

#функция расчета траекторий x, y, z, vx, vy, vz, m
def computing(x0, y0, z0, vx0, vy0, vz0, m, dt):
    R = 0.8 * m ** 0.7 * 0.0047
    x = x0 + vx0 * dt
    y = y0 + vy0 * dt
    z = z0 + vz0 * dt
    vx = np.empty_like(vx0)
    vy = np.empty_like(vy0)
    vz = np.empty_like(vz0)
    for i in range(len(m)):
        ax, ay, az = 0.0, 0.0, 0.0
        for j in range(len(m)):
            if i != j:
                dx = x[j] - x[i]
                dy = y[j] - y[i]
                dz = z[j] - z[i]
                dr = (dx * dx + dy * dy + dz * dz) ** 1.5
                if dr < (R[i] + R[j]) ** 3:
                    b = 1
                else:
                    b = 0
                ax += G * m[j] * dx / dr
                ay += G * m[j] * dy / dr
                az += G * m[j] * dz / dr
        vx[i] = vx0[i] + ax * dt
        vy[i] = vy0[i] + ay * dt
        vz[i] = vz0[i] + az * dt
    #binary center mass
    cm_x = (m[0] * x[0] + m[1] * x[1]) / (m[0] + m[1])
    cm_y = (m[0] * y[0] + m[1] * y[1]) / (m[0] + m[1])
    cm_z = (m[0] * z[0] + m[1] * z[1]) / (m[0] + m[1])
    return x, y, z, vx, vy, vz, dt, cm_x, cm_y, cm_z, b

#функция вычисления энергии звезды
def energy_computing(x, y, z, vx, vy, vz, m, cm_x, cm_y, cm_z):
    K, P = np.empty(len(m) - 1), np.empty(len(m) - 1)
    for i in range(len(m) - 1):
        V_sqrd = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]
        K[i] = m[i] * V_sqrd / 2
        dx = x[i] - cm_x
        dy = y[i] - cm_y
        dz = z[i] - cm_z
        dr = (dx * dx + dy * dy + dz * dz) ** 0.5
        P[i] = -G * m[i] / dr
    return K, P

def zero_energy(x, y, z, vx, vy, vz, m):
    cm_x = (m[0] * x[0] + m[1] * x[1] + m[2] * x[2]) / sum(m)
    cm_y = (m[0] * y[0] + m[1] * y[1] + m[2] * y[2]) / sum(m)
    cm_z = (m[0] * z[0] + m[1] * z[1] + m[2] * z[2]) / sum(m)
    K, P = np.empty(len(m)), np.empty(len(m))
    for i in range(len(m) - 1):
        V_sqrd = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]
        K[i] = m[i] * V_sqrd / 2
        for j in range(len(m)):
            if i != j:
                dx = x[i] - cm_x
                dy = y[i] - cm_y
                dz = z[i] - cm_z
                dr = (dx * dx + dy * dy + dz * dz) ** 0.5
                P[i] += -G * m[i] / dr
    K_sum, P_sum = 0, 0
    for i in range(len(m)):
        K_sum += K[i]
        P_sum += P[i]
    return K_sum + P_sum/2

#начальное состояние
x[0], y[0], z[0] = r[0] * np.cos(angle), r[0] * np.sin(angle), 0
x[1], y[1], z[1] = r[1] * np.cos(angle + np.pi), r[1] * np.sin(angle + np.pi), 0
x[2], y[2], z[2] = x_new, y_new, 0
alpha1 = np.arctan((y[0] - y[2]) / (x[0] - x[2]))
alpha2 = np.arctan((y[1] - y[2]) / (x[1] - x[2]))
acc1 = G * m[2] / ((y[0] - y[2]) ** 2 + (x[0] - x[2]) ** 2)
acc2 = G * m[2] / ((y[1] - y[2]) ** 2 + (x[1] - x[2]) ** 2)
vx[0], vy[0], vz[0] = np.sqrt(G * (sum(m) - m[2]) / (r[0] + r[1])) / (1 + m[0]/m[1]) * np.cos(angle + np.pi/2) - acc1 * dt * np.cos(alpha1), np.sqrt(G * (sum(m) - m[2]) / (r[0] + r[1])) / (1 + m[0]/m[1]) * np.sin(angle + np.pi/2) - acc1 * dt * np.sin(alpha1), 0
vx[1], vy[1], vz[1] = np.sqrt(G * (sum(m) - m[2]) / (r[0] + r[1])) / (1 + m[1]/m[0]) * np.cos(angle + 3/2*np.pi) - acc2 * dt * np.cos(alpha2), np.sqrt(G * (sum(m) - m[2]) / (r[0] + r[1])) / (1 + m[1]/m[0]) * np.sin(angle + 3/2*np.pi) - acc2 * dt * np.sin(alpha2), 0
vx[2], vy[2], vz[2] = -v, 0, 0
cm_x, cm_y, cm_z = 0, 0, 0

#рабочий цикл
if zero_energy(x, y, z, vx, vy, vz, m) >= 0:
    flag = 0
    i = 0
    while i < steps:
        K, P = energy_computing(x, y, z, vx, vy, vz, m, cm_x, cm_y, cm_z)
        x, y, z, vx, vy, vz, dt, cm_x, cm_y, cm_z, b = computing(x, y, z, vx, vy, vz, m, dt)
        #flyby
        if ((x[2] - cm_x) ** 2 + (y[2] - cm_y) ** 2 + (z[2] - cm_z) ** 2) >= (x_new ** 2 + y_new ** 2) and flag == 0 and i > 5000:
            print('flyby')
            flag = 1
            steps = i + 2500
        #close objects
        if b == 1:
            print('objects were too close')
            shutil.os.remove(f'{filename}.npz')
            break
        #runaway or exchanged
        if (K[0] + P[0] + K[1] + P[1]) > 0:
            if (vx[0] ** 2 + vy[0] ** 2 + vz[0] ** 2 >= (18.6/1731.48148) ** 2 and (flag == 0)) or (vx[1] ** 2 + vy[1] ** 2 + vz[1] ** 2 >= (18.6/1731.48148) ** 2 and (flag == 0)):
                print('binary becomes runaway')
                steps = i + 6000
                flag = 1
            elif (vx[0] ** 2 + vy[0] ** 2 + vz[0] ** 2 < (18.6/1731.48148) ** 2 and (flag == 0)) or (vx[1] ** 2 + vy[1] ** 2 + vz[1] ** 2 < (18.6/1731.48148) ** 2 and (flag == 0)):
                print('binary exchanged')
                steps = i + 6000
                flag = 1
        #saving
        if (i + 1) == 100:
            save(filename, x, y, z, vx, vy, vz, m, round(dt * (i + 1) / 365.2522, 2), K, P)
        elif (i + 1) % 100 == 0 and not((i + 1) - 100 == 0):
            if bool(os.path.exists(f'{filename}.npz')) == 1:
                imported = np.load(f'{filename}.npz')
                imported = dict(imported)
                imported['x'] = np.vstack([imported['x'], x])
                imported['y'] = np.vstack([imported['y'], y]) 
                imported['z'] = np.vstack([imported['z'], z])
                imported['vx'] = np.vstack([imported['vx'], vx])
                imported['vy'] = np.vstack([imported['vy'], vy])
                imported['vz'] = np.vstack([imported['vz'], vz])
                imported['t'] = np.append(imported['t'], np.array(round(dt * (i + 1) / 365.2522, 2)))
                imported['P'] = np.vstack([imported['P'], P])
                imported['K'] = np.vstack([imported['K'], K])
                save(filename, imported['x'], imported['y'], imported['z'], imported['vx'], imported['vy'], imported['vz'], imported['m'], imported['t'], imported['K'], imported['P'])
            else:
                print("File doesn't exist")
                break
        i += 1
else:
    print('System will be stable')
