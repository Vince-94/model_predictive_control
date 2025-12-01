#!/usr/bin/env python3

import pytest
from test.env_fixtures import sys_model_fixture, mpc_params_fixture

from system_discretization import SystemDiscretization
from mpc_dataclass import SystemModelData

import numpy as np


class TestSystemDiscretization():

    def test_euler_discretization(self, sys_model_fixture, mpc_params_fixture):
        self.Ad = np.array([
            [1, 0, 0.1, 0],
            [0, 1, 0, 0.1],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])
        self.Bd = np.array([
            [0.1, 0],
            [0, 0.1],
            [0.1, 0],
            [0, 0.1]])
        self.Cd = np.array([
            [0.1, 0, 0, 0],
            [0, 0.1, 0, 0],
            [0, 0, 0.1, 0]])
        self.Dd = np.array([
            [0, 0],
            [0, 0],
            [0, 0]])

        disc_system = SystemDiscretization(cont_sys_model=sys_model_fixture, Ts=mpc_params_fixture.Ts).euler_discretization()

        np.testing.assert_allclose(disc_system.A, self.Ad, verbose=True)
        np.testing.assert_allclose(disc_system.B, self.Bd, verbose=True)
        np.testing.assert_allclose(disc_system.C, self.Cd, verbose=True)
        np.testing.assert_allclose(disc_system.D, self.Dd, verbose=True)

