import os
import matplotlib.pyplot as plt
import numpy as np

a = np.logspace(-1, 1, 10)
def file_reading(M, a):
    os.chdir(f'C:\\Science\\{M}M')
    data_dict = {}
    with open('output.txt', 'r') as f:
        data = f.readlines()
        for i in range(len(data)):
            data[i] = data[i][:-1]
            try:
                if 'system' in data[i] and 'system' not in data[i + 1] and 'objects' not in data[i+1] and 'objects' not in data[i+2]:
                    data_dict[data[i]] = data[i + 1]
            except IndexError:
                data_dict[data[i]] = data[i + 1]
                print('error')
    results = [[], [], [], [], [], [], [], [], [], []]
    for i in data_dict:
        number = (int(i.replace('system_state', '')) - 1) // 100
        results[number].append(data_dict[i][:-1])
    fl_part, ra_part, ex_part = [], [], []
    for i in range(len(results)):
        fl_part.append(results[i].count('flyby')/len(results[i]))
        ra_part.append(results[i].count('binary becomes runaway')/len(results[i]))
        ex_part.append(results[i].count('binary exchanged')/len(results[i]))
    plt.figure()
    plt.xlabel('a [AU]')
    plt.ylabel('Branching ratio')
    if M == 1:
        plt.title(r'Binary, 1 $M_{sun}$')
    if M == 1.4:
        plt.title(r'Binary, 1.4 $M_{sun}$')
    if M == 6:
        plt.title(r'Binary, 6 $M_{sun}$')
    if M == 10:
        plt.title(r'Binary, 10 $M_{sun}$')
    plt.plot(a, fl_part, label = 'flybys')
    plt.plot(a, ra_part, label = 'Runaways')
    plt.plot(a, ex_part, label = 'Exchanges')
    plt.legend(loc = 'best')
    os.chdir('C:\\Users\\nmago\\Desktop')
    plt.savefig(f'{M}_masses.pdf')
    plt.show()
    return
file_reading(1, a)
file_reading(1.4, a)
file_reading(6, a)
file_reading(10, a)
