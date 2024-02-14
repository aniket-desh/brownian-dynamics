import numpy as np
import matplotlib.pyplot as plt

# get the 4th line of input.txt
def get_nbb(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        return int(lines[4].split()[-1])

# get the 1st column of stats.txt and save as a numpy array
def get_time(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        return np.array([float(line.split()[0]) for line in lines])

# get the 4th column of stats.txt and save as a numpy array
def get_energy(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        return np.array([float(line.split()[3]) for line in lines])

number_backbone = get_nbb('input.txt')
time = get_time('stats.txt')
energy = get_energy('stats.txt')

print('Number of backbone beads:', number_backbone)
print('Time:', time)
print('Energy:', energy)
plt.plot(time, energy / number_backbone)
plt.show()

