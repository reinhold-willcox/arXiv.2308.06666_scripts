import numpy as np

mm1 = np.append(np.linspace(2, 100, 99), 6.3)


with open('grid.txt', 'w') as f:
    for m1 in mm1:
        f.write("--initial-mass-1 {} --initial-mass-2 0.5 -a 1000\n".format(m1))



