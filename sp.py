import sys
import numpy as np
from numpy import linalg as la


# Atomic masses and van der Waals radii
num_list = {'1 ': 'H ', '6 ': 'C ', '7 ': 'N ', '8 ': 'O ', '9 ': 'F ', '14': 'Si',
            '15': 'P ', '16': 'S ', '17': 'Cl', '35': 'Br', '53': 'I '}
m_list = {'H ': 1.008, 'C ': 12.010, 'N ': 14.006, 'O ': 15.999, 'F ': 18.998, 'Si': 28.086,
          'P ': 30.974, 'S ': 32.059, 'Cl': 35.446, 'Br': 79.901, 'I ': 126.904}
r_list = {'H ': 1.20, 'C ': 1.77, 'N ': 1.66, 'O ': 1.50, 'F ': 1.46, 'Si': 2.19,
          'P ': 1.90, 'S ': 1.89, 'Cl': 1.82, 'Br': 1.86, 'I ': 2.04}

# Reading the input coordinates of atoms
File = sys.argv[1]
if File[-3:] == 'sdf':
    with open(File) as file:
        Lines = file.readlines()
    num_atoms = int(Lines[3][0:3])
    coords = np.ndarray((num_atoms, 3), dtype=float)
    masses = np.array((num_atoms, 1), dtype=float)
    elements = np.array([Lines[4 + a][31:33] for a in range(num_atoms)])
    for a in range(num_atoms):
        coords[a, 0] = Lines[4 + a][0:10]
        coords[a, 1] = Lines[4 + a][11:20]
        coords[a, 2] = Lines[4 + a][21:30]
elif File[-3:] == 'xyz':
    with open(File) as file:
        Lines = file.readlines()
    num_atoms = len(Lines)
    coords = np.ndarray((num_atoms, 3), dtype=float)
    masses = np.array((num_atoms, 1), dtype=float)
    if Lines[0][0].isdigit():
        elements = np.array([num_list[Lines[a][0:2]] for a in range(num_atoms)])
    else:
        elements = np.array([Lines[a][0:2] for a in range(num_atoms)])
    for a in range(num_atoms):
        coords[a, 0] = Lines[a][8:20]
        coords[a, 1] = Lines[a][25:37]
        coords[a, 2] = Lines[a][42:54]
else:
    print('Unsupported format. Please provide coordinates as .sdf or .xyz file.')
    sys.exit()
masses = np.array([m_list[elements[a]] for a in range(num_atoms)])

# Determination of the coordinates of the center of mass of a molecule
mol_mass = np.sum(masses)
com = np.array([np.sum([coords[a, axis] * masses[a] for a in range(num_atoms)]) / mol_mass for axis in range(3)])

# Translation of the coordinates of all atoms to the center of mass
for a in range(num_atoms):
    coords[a, :] = coords[a, :] - com

# Determination of the moment of inertia tensor, I
x = coords[:, 0]
y = coords[:, 1]
z = coords[:, 2]
I = np.ndarray((3, 3))
I[0, 0] = np.sum(masses * (y ** 2 + z ** 2))
I[0, 1] = -np.sum(masses * x * y)
I[0, 2] = -np.sum(masses * x * z)
I[1, 0] = I[0, 1]
I[1, 1] = np.sum(masses * (x ** 2 + z ** 2))
I[1, 2] = -np.sum(masses * y * z)
I[2, 0] = I[0, 2]
I[2, 1] = I[1, 2]
I[2, 2] = np.sum(masses * (x ** 2 + y ** 2))

# Calculation of the eigenvectors of the inertia matrix, i.e. and principal rotation axes of the molecule
eigenvalues, eigenvectors = la.eig(I)

# Using the rotation matrix to orient the axes of the coordinate system along the principal axes of the molecule
coords = np.matmul(coords, eigenvectors)

# Drawing a sphere around each atom with a radius equal to its van der Waals radius
spheres = np.empty((0, 3), dtype=float)
for a in range(num_atoms):
    radius = r_list[elements[a]]
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    for t in u:
        for p in v:
            x = coords[a, 0] + radius * np.cos(t) * np.sin(p)
            y = coords[a, 1] + radius * np.sin(t) * np.sin(p)
            z = coords[a, 2] + radius * np.cos(p)
            sphere = np.array([x, y, z]).transpose()
            spheres = np.vstack((spheres, sphere))

# Calculation of geometric parameters of the molecule: length, width and thickness
params = [np.max(spheres[:, 0]) - np.min(spheres[:, 0]), np.max(spheres[:, 1]) - np.min(spheres[:, 1]),
          np.max(spheres[:, 2]) - np.min(spheres[:, 2])]
params.sort()
length = params[2]
width = params[1]
thickness = params[0]

# Calculation of the sphericity parameter by dividing the thickness of the molecule by its length
sp = thickness / length
print(File[:-4], 'sp =', np.round(sp, 2))
