#!/usr/bin/env python3

import pytest
from test.env_fixtures import sys_model_fixture, mpc_params_fixture, signals_fixture

from lti_mpc import LinearTimeIvariantMPC


class TestLinearTimeIvariantMPC():

    def test_on_step(self, sys_model_fixture, mpc_params_fixture, signals_fixture):
        lti_mpc = LinearTimeIvariantMPC(sys_model_fixture, mpc_params_fixture)
        u_star = lti_mpc.on_step(signals_fixture)
        assert len(u_star) == mpc_params_fixture.nu
