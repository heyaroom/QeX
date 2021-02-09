# caution! this method is authorized (officially cannot be used)

import numpy as np

def optimize(model, p_seed, iteration, initp=None):
    m = model
    n_param = m.n_param
    gs_min_width_relative = False
    gs_min_width = 0.001
    gs_similarity_threshold = 0.5

    np.random.seed(p_seed)
    if initp is None:
        param = np.random.random(size=n_param) * 2 * np.pi
    else:
        param = initp
        
    recycle_z0 = None
    gs_tool = {'param': None, 'z0': np.inf, 'shift': None}  # For Golden Section Search

    m.step(np.copy(param), evaluate=True)

    for i in range(iteration):

        idx = i % m.n_param

        if recycle_z0 is not None:
            z0 = recycle_z0
            recycle_z0 = None
        else:
            z0 = m.step(phi=np.copy(param))

        if idx == 0:
            gs_tool['param'] = np.copy(param)
            gs_tool['z0'] = z0

        p = np.copy(param)
        p[idx] = param[idx] + np.pi / 2
        z1 = m.step(phi=p)

        p = np.copy(param)
        p[idx] = param[idx] - np.pi / 2
        z3 = m.step(phi=p)

        z2 = z1 + z3 - z0

        c = (z0 + z1 + z2 + z3) / 4
        a = np.sqrt((z0 - z2) ** 2 + (z1 - z3) ** 2) / 2
        b = np.arctan((z1 - z3) / ((z0 - z2) + 1e-32 * (z0 == z2))) + param[idx]
        b += 0.5 * np.pi + 0.5 * np.pi * np.sign((z0 - z2) + 1e-32 * (z0 == z2))

        param[idx] = b
        recycle_z0 = c - a
        
        if idx == m.n_param - 1:

            shift = param - gs_tool['param']

            if gs_tool['shift'] is not None:

                eps = 1e-12

                norm_shift = np.linalg.norm(shift)
                norm_shift += (norm_shift == 0) * eps  # Avoid Zero Division

                norm_old_shift = np.linalg.norm(gs_tool['shift'])
                norm_old_shift += (norm_old_shift == 0) * eps  # Avoid Zero Division

                #similarity = (shift / norm_shift) @ (gs_tool['shift'] / norm_old_shift)
                similarity = np.dot((shift / norm_shift),(gs_tool['shift'] / norm_old_shift))

                if similarity > gs_similarity_threshold:
                    if gs_min_width_relative:
                        tolerance = norm_shift * gs_min_width
                    else:
                        tolerance = gs_min_width

                    param, recycle_z0 = golden_section_search(
                        m.step,
                        gs_tool['param'], gs_tool['z0'],
                        param, recycle_z0,
                        tolerance=tolerance
                    )
                    shift = None

            gs_tool['shift'] = shift

        try:
            m.step(np.copy(param), evaluate=True, stop_signal=True)
        except:
            break
def golden_section_search(func, param1, loss1, param2, loss2=None, tolerance=1e-10, max_iter=100):

    if loss1 is None:
        loss1 = func(param1)

    if loss2 is None:
        loss2 = func(param2)

    # Set  loss_left > loss_mid
    if loss1 < loss2:
        loss_left, loss_mid = loss2, loss1
        param_left, param_mid = param2, param1
    else:
        loss_left, loss_mid = loss1, loss2
        param_left, param_mid = param1, param2

    shift = param_mid - param_left

    if np.linalg.norm(shift) == 0:
        return param_mid, loss_mid

    golden = (1 + np.sqrt(5)) / 2

    # Expansion Part
    while True:
        param_right = param_mid + shift * golden
        loss_right = func(param_right)

        if loss_right >= loss_mid:
            break

        param_left, param_mid = param_mid, param_right
        loss_left, loss_mid = loss_mid, loss_right
        shift *= golden

    # Partition Part
    big_right = True
    while True:
        if np.linalg.norm(shift * (1 + golden)) < tolerance:  # width < tolerance
            break

        param_proposal = param_mid + shift * golden if big_right else param_mid - shift * golden
        loss_proposal = func(param_proposal)

        if loss_proposal < loss_mid:
            param_mid = param_proposal
            loss_mid = loss_proposal
        else:
            big_right = not big_right

        shift /= golden

    return param_mid, loss_mid

