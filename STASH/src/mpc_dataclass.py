#!/usr/bin/env python3
'''

'''
from dataclasses import dataclass
import numpy as np

@dataclass
class SignalsData:
    x: np.array = None  # size(nx)
    ref: np.array = None # size(ny)
    u: np.array = None # size(ny)
    x_dot: np.array = None  # size(nx)
    y: np.array = None # size(ny)


@dataclass
class SystemModelData:
    A: np.array # size(nx, nx)
    B: np.array # size(nx, nu)
    C: np.array # size(ny, nu)
    D: np.array # size(ny, nu)


@dataclass
class PredictionModelData:
    A_bar: np.array # size(nx Np, nx)
    B_bar: np.array # size(nx Np, nu Nc)
    C_bar: np.array # size(ny Np, nx)
    D_bar: np.array # size(ny Np, nu Nc)
    REF: np.array # size(ny Np)
    U0: np.array # size(nu Nc)
    Q_bar: np.array # size(ny Np, ny Np)
    R_bar: np.array # size(nu Np, nu Np)
    X_c_min: np.array # size(nx Np, 1)
    X_c_max: np.array # size(nx Np, 1)
    U_c_min: np.array # size(nu Nc, 1)
    U_c_max: np.array # size(nu Nc, 1)
    dU_c_min: np.array # size(nu Nc, 1)
    dU_c_max: np.array # size(nu Nc, 1)


@dataclass
class QuadraticProblemData:
    H: np.array # size(nu Np, nu Np)
    f: np.array # size(nu Np, ?)
    A_ineq: np.array # size()
    b_ineq: np.array # size()
    A_eq: np.array = None # size()
    b_eq: np.array = None # size()


@dataclass
class MpcParamData:
    nx: int # size(A, 0)
    ny: int # size(C, 0)
    nu: int # size(B, 1)
    Ts: float
    Np: int
    Nc: int
    Q: np.array # size(nx,nx)
    R: np.array # size(nu,nu)
    x_c_min: np.array # size(nx,1)
    x_c_max: np.array # size(nx,1)
    u_c_min: np.array # size(nu,1)
    u_c_max: np.array # size(nu,1)
    du_c_min: np.array # size(nu,1)
    du_c_max: np.array # size(nu,1)


if __name__ == '__main__':
    print('MPC Dataclass')
