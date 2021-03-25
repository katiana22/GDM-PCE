import numpy as np


def function(x, y, r0=0.5, eps_in=3, eps_out=1, efield=10):
    """
    Dielectric cylinder of radius 'r0' is suspended in a homogeneous electric field
    with strength 'efield'. The permittivity of the cylinder's material is
    'eps_in'. The permittivity of the surrounding space is 'eps_out'.
    The electric potential is given by 'analytical_potential'

    """
    radius = np.sqrt(x ** 2 + y ** 2)
    eps_ratio = eps_in / eps_out
    k = ((eps_ratio - 1) / (eps_ratio + 1)) * r0 ** 2
    if radius > r0:
        return -efield * x * (1 - k / (radius ** 2))
    else:
        return -efield * x * 2 / ((eps_in / eps_out) + 1)

