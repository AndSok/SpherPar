# SpherPar
SpherPar is a Python script for the calculation of the sphericity parameter according to the procedure described in []().

# Usage
Provide a file containing the geometry of a molecule in .sdf or .xyz format:
```
python sp.py molecule.sdf
```
or
```
python sp.py molecule.xyz
```
The output is a sphericity parameter of the desired molecule:
```
python sp.py benzene.sdf
benzene sp = 0.48
```
```
python sp.py n-eicosane.xyz
n-eicosane sp = 0.15
```

# How to cite
Please cite the following papers when using SpherPar:\
1.
