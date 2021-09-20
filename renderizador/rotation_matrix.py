import numpy as np
from math import sin, cos, pi

def sqr(a, b):
    return np.round(1 - (2*((a**2) + (b**2))), 10)

def neg(a, b, c, d):
    return np.round(2*((a*b) - (c*d)), 10)

def pos(a, b, c, d):
    return np.round(2*((a*b) + (c*d)), 10)

def rotation_matrix(xyzt):

    #quatern"i"on = [cos(xyzt[3]/2), s"i"n(xyzt[3]/2)*xyzt[0], s"i"n(xyzt[3]/2)*xyzt[1], s"i"n(xyzt[3]/2)*xyzt[2]]

    q = {
        "r":cos(xyzt[3]/2),
        "i":sin(xyzt[3]/2)*xyzt[0],
        "j":sin(xyzt[3]/2)*xyzt[1],
        "k":sin(xyzt[3]/2)*xyzt[2],
    }

    return np.array([
            [sqr(q["j"], q["k"]),                   neg(q["i"], q["j"], q["k"], q["r"]),    pos(q["i"], q["k"], q["j"], q["r"]),    0],
            [pos(q["i"], q["j"], q["k"], q["r"]),   sqr(q["i"], q["k"]),                    neg(q["j"], q["k"], q["i"], q["r"]),    0],
            [neg(q["i"], q["k"], q["j"], q["r"]),   pos(q["j"], q["k"], q["i"], q["r"]),    sqr(q["i"], q["j"]),                    0],
            [0,                                     0,                                      0,                                      1],
        ])
    
# vec = [1.0, 0.0, 0.0, -1.57]

# matriz = rotation_matrix(vec)

# print(matriz)