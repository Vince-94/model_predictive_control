#!/usr/bin/env python3

import pytest

from mpc_dataclass import SystemModelData, MpcParamData, SignalsData

import numpy as np


@pytest.fixture
def sys_model_fixture():
    A = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0]])
    B = np.array([
        [1, 0],
        [0, 1],
        [1, 0],
        [0, 1]])
    C = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0]])
    D = np.array([
        [0, 0],
        [0, 0],
        [0, 0]])
    return SystemModelData(A=A, B=B, C=C, D=D)


@pytest.fixture
def mpc_params_fixture():
    nx = 4
    ny = 3
    nu = 2
    Ts = 0.1
    Np = 3
    Nc = 2
    Q = np.diag([1, 1, 1])
    R = np.diag([1, 1])
    x_c_min = np.array([
        [-np.inf],
        [-np.inf],
        [-np.inf],
        [-np.inf]])
    x_c_max = np.array([
        [np.inf],
        [np.inf],
        [np.inf],
        [np.inf]])
    u_c_min = np.array([
        [-np.inf],
        [-np.inf]])
    u_c_max = np.array([
        [np.inf],
        [np.inf]])
    du_c_min = np.array([
        [-np.inf],
        [-np.inf]])
    du_c_max = np.array([
        [np.inf],
        [np.inf]])
    return MpcParamData(
        nx=nx, ny=ny, nu=nu, Ts=Ts, Np=Np, Nc=Nc, Q=Q, R=R,
        x_c_min=x_c_min, x_c_max=x_c_max,
        u_c_min=u_c_min, u_c_max=u_c_max,
        du_c_min=du_c_min, du_c_max=du_c_max)


@pytest.fixture
def signals_fixture():
    x0 = np.array([[0], [0], [0], [0]])
    ref = np.array([[0], [0], [0]])
    u0 = np.array([[0], [0]])
    return SignalsData(x=x0, ref=ref, u=u0)
