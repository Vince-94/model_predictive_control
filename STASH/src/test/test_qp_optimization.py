#!/usr/bin/env python3

import pytest
from test.env_fixtures import sys_model_fixture, mpc_params_fixture, signals_fixture

from system_prediction import SystemPrediction
from qp_optimization import QuadraticProgramOptimization
from mpc_dataclass import SystemModelData, MpcParamData, SignalsData

import numpy as np


class TestQuadraticProgramOptimization():

    @pytest.fixture
    def sys_pred_fixture(self, sys_model_fixture, mpc_params_fixture, signals_fixture):
        sys_pred =  SystemPrediction(sys_model=sys_model_fixture, mpc_params=mpc_params_fixture)
        sys_pred.prediction_invariant_matrices()
        sys_pred.prediction_varying_matrices(signals_fixture)
        return sys_pred.get_prediction_matrices


    def test_qp_matrices(self, mpc_params_fixture, signals_fixture, sys_pred_fixture):
        self.qp_opt = QuadraticProgramOptimization(mpc_params_fixture.nx, mpc_params_fixture.nu, mpc_params_fixture.Np, mpc_params_fixture.Nc, sys_pred_fixture)
        qp_matr = self.qp_opt.qp_matrices(signals_fixture.x, sys_pred_fixture)

        assert qp_matr.H.shape == (mpc_params_fixture.nu*mpc_params_fixture.Nc, mpc_params_fixture.nu*mpc_params_fixture.Nc)
        assert qp_matr.f.shape == (1, mpc_params_fixture.nu*mpc_params_fixture.Nc)
        assert qp_matr.A_ineq.shape == (2*mpc_params_fixture.nx*mpc_params_fixture.Np+4*mpc_params_fixture.nu*mpc_params_fixture.Nc, mpc_params_fixture.nu*mpc_params_fixture.Nc)
        assert qp_matr.b_ineq.shape == (2*mpc_params_fixture.nx*mpc_params_fixture.Np+4*mpc_params_fixture.nu*mpc_params_fixture.Nc, 1)

        assert np.allclose(qp_matr.H, qp_matr.H.T), "Matrix H is not symmetric"


    def test_qp_solve(self, mpc_params_fixture, signals_fixture, sys_pred_fixture):
        self.qp_opt = QuadraticProgramOptimization(mpc_params_fixture.nx, mpc_params_fixture.nu, mpc_params_fixture.Np, mpc_params_fixture.Nc, sys_pred_fixture)
        self.qp_matr = self.qp_opt.qp_matrices(signals_fixture.x, sys_pred_fixture)
        U_star = self.qp_opt.qp_solve(signals_fixture.x, sys_pred_fixture)

        assert len(U_star) == mpc_params_fixture.nu*mpc_params_fixture.Nc


    #todo this could be in another class
    def test_qp_solver_quadprog(self, mpc_params_fixture, signals_fixture, sys_pred_fixture):
        self.qp_opt = QuadraticProgramOptimization(mpc_params_fixture.nx, mpc_params_fixture.nu, mpc_params_fixture.Np, mpc_params_fixture.Nc, sys_pred_fixture)
        qp_matr = self.qp_opt.qp_matrices(signals_fixture.x, sys_pred_fixture)
        U_star = self.qp_opt.qp_solver_quadprog(qp_matr.H, qp_matr.f, qp_matr.A_ineq, qp_matr.b_ineq)

        # np.testing.assert_array_almost_equal(xf, [0.4761905, 1.0476190, 2.0952381])
        # np.testing.assert_almost_equal(f, -2.380952380952381)
        # np.testing.assert_almost_equal(xu, [0, 5, 0])
        # np.testing.assert_array_equal(iters, [3, 0])
        # np.testing.assert_array_almost_equal(lagr, [0.0000000, 0.2380952, 2.0952381])

