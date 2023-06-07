import os
import numpy as np

path = "C:\\Science\\10M"
os.chdir(path)

with open ('parallel.txt', 'w') as f:
    koef = 1731.48148
    dt = np.logspace(-2.5, -0.3, 10)
    a = np.logspace(-1, 1, 10)
    name = 0
    f.write('#!/usr/bin/env python\n\n')
    for j in range(len(a)):
        for i in range(100):
            y = np.sqrt(np.random.uniform(0, 25 * a[j] ** 2))
            v = abs(np.random.normal(0, 5.5) / koef)
            name += 1
            f.write(f'#a = {a[j]}, m = 10, x = {30 * a[j]}, y = {y}, v = {v}, steps = 250000, dt = {dt[j]}, filename = system_state{name}\n')
            f.write(f'/usr/bin/python3 ../modelling_with_energy.py {a[j]} 10 {30 * a[j]} {y} {v} 250000 {dt[j]} system_state{name}\n')
            f.write('\n')
