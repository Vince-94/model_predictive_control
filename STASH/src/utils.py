#!/usr/bin/env python3

from mpc_dataclass import SystemModelData, MpcParamData, SignalsData

import numpy as np

def check_system_dimensions(system_model: SystemModelData, mpc_params: MpcParamData, signals: SignalsData) -> bool:
    '''
    Check if MPC parameters (nx, ny, nu) match the system model dimensions
    '''
    assert mpc_params.nx == np.size(system_model.A, 0), "nx set wrong!"
    assert mpc_params.ny == np.size(system_model.C, 0), "ny set wrong!"
    assert mpc_params.nu == np.size(system_model.B, 1), "nu set wrong!"

    block_up = np.hstack((
        system_model.A, system_model.B
    ))
    block_down = np.hstack((
        system_model.C, system_model.D
    ))
    block = np.vstack((
        block_up,
        block_down
    ))
    assert block.shape == (mpc_params.nx + mpc_params.ny, mpc_params.nx + mpc_params.nu)

    assert mpc_params.nx == len(signals.x), f"size(nx)={mpc_params.nx} != size(x)={len(signals.x)}"
    assert mpc_params.ny == len(signals.ref), f"size(ny)={mpc_params.ny} != size(ref)={len(signals.ref)}"
    assert mpc_params.nu == len(signals.u), f"size(nu)={mpc_params.nu} != size(u)={len(signals.u)}"

    assert mpc_params.Q.shape == (mpc_params.ny, mpc_params.ny)
    assert mpc_params.R.shape == (mpc_params.nu, mpc_params.nu)

    assert len(mpc_params.x_c_min) == mpc_params.nx
    assert len(mpc_params.x_c_max) == mpc_params.nx
    assert len(mpc_params.u_c_min) == mpc_params.nu
    assert len(mpc_params.u_c_min) == mpc_params.nu
    assert len(mpc_params.du_c_max) == mpc_params.nu
    assert len(mpc_params.du_c_max) == mpc_params.nu


if __name__ == '__main__':
    print('Utility functions')
