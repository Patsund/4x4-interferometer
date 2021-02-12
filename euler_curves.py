import numpy as np
from scipy.special import fresnel
from scipy.optimize import brentq

from gdshelpers.parts.waveguide import Waveguide


def rotate(point, angle, origin):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point
    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy


def EulerBendPoints_Simple(angle_amount=np.pi / 2., radius=10.0, clockwise=False, resolution=200.0):
    """ Euler bend, no transformation, emerging from the origin. Given as coordinate points"""

    # End angle

    # clockwise = bool((angle_amount % (2 * np.pi)) > np.pi)

    if clockwise:
        eth = (-angle_amount) % (2 * np.pi)
    else:
        eth = angle_amount % (2 * np.pi)

    # If bend is trivial, return a trivial shape

    if np.isclose(eth, 0.0):
        return [np.array((0, 0))]

    # Curve min radius
    R = radius

    # Total displaced angle
    th = eth / 2.0

    # Total length of curve
    Ltot = 4 * R * th

    ###################
    ## Compute curve ##
    ###################

    a = np.sqrt(R ** 2. * th)
    sq2pi = np.sqrt(2. * np.pi)

    # Function for computing curve coords
    (fasin, facos) = fresnel(np.sqrt(2. / np.pi) * R * th / a)

    def XY(s):
        if th == 0:
            return np.array((0.0, 0.0))
        elif s <= Ltot / 2:
            (fsin, fcos) = fresnel(s / (sq2pi * a))
            X = sq2pi * a * fcos
            Y = sq2pi * a * fsin
        else:
            (fsin, fcos) = fresnel((Ltot - s) / (sq2pi * a))
            X = sq2pi * a * (facos
                             + np.cos(2 * th) * (facos - fcos)
                             + np.sin(2 * th) * (fasin - fsin))
            Y = sq2pi * a * (fasin
                             - np.cos(2 * th) * (fasin - fsin)
                             + np.sin(2 * th) * (facos - fcos))

        return np.array((X, Y))

    # Parametric step size
    # print('th, res:', th, resolution)
    step = Ltot / (int(th * resolution) + 1)

    # Generate points
    points = []
    for i in range(0, int(round(Ltot / step)) + 1):
        pt = XY(i * step)
        if clockwise:
            points.append(np.array([pt[0], - pt[1]]))
        else:
            points.append(pt)

    return points


def EulerBendPoints(start_point=(0, 0),
                    radius=10.,
                    input_angle=0.0,
                    output_angle=np.pi / 2.,
                    clockwise=False,
                    resolution=200.0):
    """ Euler bend with a given start point, input angle, output angle and radius"""
    pts0 = EulerBendPoints_Simple(angle_amount=output_angle - input_angle, radius=radius, clockwise=clockwise,
                                  resolution=resolution)
    origin_pt = np.array(start_point)
    pts = [rotate(pt, input_angle, (0, 0)) + origin_pt for pt in pts0]
    return pts


def EulerBendPoints_Relative(start_point=(0, 0),
                             radius=10.,
                             input_angle=0.0,
                             angle_amount=np.pi / 2.,
                             clockwise=False,
                             resolution=200.0):
    """Relative Euler bend: bend with relative turning angle instead of absolute end angle. JWS Sept 2015"""
    return EulerBendPoints(start_point=start_point, radius=radius, input_angle=input_angle,
                           output_angle=input_angle + angle_amount, clockwise=clockwise, resolution=resolution)


def EulerEndPt(start_point=(0.0, 0.0),
               radius=10.0,
               input_angle=0.0,
               clockwise=False,
               angle_amount=np.pi / 2.):
    """Gives the end point of a simple Euler bend"""

    if clockwise:
        eth = (-angle_amount) % (2 * np.pi)
    else:
        eth = angle_amount % (2 * np.pi)

    th = eth / 2.0
    R = radius
    # clockwise = bool(angle_amount < 0)

    (fsin, fcos) = fresnel(np.sqrt(2 * th / np.pi))

    a = 2 * np.sqrt(2 * np.pi * th) * (np.cos(th) * fcos + np.sin(th) * fsin)
    r = a * R
    X = r * np.cos(th)
    Y = r * np.sin(th)

    if clockwise:
        Y *= -1

    pt = np.array((X, Y)) + np.array(start_point)
    pt = rotate(point=pt, angle=input_angle, origin=start_point)

    return pt


def EulerLength(radius=10.0,
                angle_amount=np.pi / 2.,
                clockwise=False):
    if clockwise:
        eth = (-angle_amount) % (2 * np.pi)
    else:
        eth = angle_amount % (2 * np.pi)
    th = eth / 2.0
    return 4 * radius * th


def EulerSBend_Points(start_point=(0.0, 0.0),
                      offset=5.0,
                      radius=10.0,
                      input_angle=0.0,
                      resolution=200.):
    """An Euler s-bend with parallel input and output, separated by an offset."""

    init_point = np.array(start_point)

    # Check that a shape should be drawn
    if offset == 0:
        return [(0, 0)]

    in_clockwise = bool(offset < 0)

    # Function to find root of
    def froot(th):
        return 2 * EulerEndPt((0., 0.), radius, 0., False, th)[1] - abs(offset)

    # Get direction
    if offset >= 0:
        dir = +1
    else:
        dir = -1

    # Check whether offset requires straight section
    a = 0.0
    b = np.pi / 2.
    fa = froot(a)
    fb = froot(b)

    if fa * fb < 0:
        # Offset can be produced just by bends alone
        angle = dir * brentq(froot, 0., np.pi / 2.)
        extra_y = 0.0
    else:
        # Offset is greater than max height of bends
        angle = dir * np.pi / 2.
        extra_y = - dir * fb

    # First bend
    points = EulerBendPoints_Relative(init_point, radius, 0, angle, in_clockwise, resolution)
    # Second bend
    points = np.concatenate(
        (
            points, EulerBendPoints_Relative(points[-1] + np.array((0, extra_y)), radius, angle, -angle,
                                             (not in_clockwise), resolution)))

    # Check
    if abs(points[-1][1] - offset - start_point[1]) > 0.00005:
        print("WARNING: EulerSBend: Output curve has incorrect offset by %f um" % (points[-1][1] - offset))

    points = [rotate(pt, input_angle, init_point) for pt in points]
    return points


def EulerWiggle_Points(start_point=(0.0, 0.0),
                       input_angle=0.0,
                       radius=10.0,
                       target_path_length=None,
                       target_crow_length=None,
                       internal_angle_mod=0.0,
                       N_turns=10,
                       mirrored=False,
                       resolution=200,
                       inout_space=1):
    """Euler-bend-based wiggle, with N_turns lobes, each subtending 180+internal_angle_mod degrees."""

    # Check inputs
    if type(N_turns) is not int or N_turns < 1:
        raise AttributeError('Parameter N_turns must be of type int, and greater than 0.')

    ################
    # CALCULATIONS #
    ################

    def get_lens(radius, internal_angle_mod, N_turns):
        io_bend = EulerBendPoints_Relative(
            radius=radius,
            angle_amount=np.pi / 2. + internal_angle_mod,
            clockwise=False)

        wiggle_bend = EulerBendPoints_Relative(
            radius=radius,
            input_angle=(np.pi / 2. + internal_angle_mod),
            angle_amount=-(np.pi + 2 * internal_angle_mod),
            clockwise=True)

        # Get crow length (start-to-finish as the crow flies)
        cl_io = io_bend[-1][0]
        cl_wiggle = wiggle_bend[-1][0]

        # Get path length
        l_io = EulerLength(radius=radius, angle_amount=(np.pi / 2. + internal_angle_mod), clockwise=False)
        l_wiggle = EulerLength(radius=radius, angle_amount=-(np.pi + 2 * internal_angle_mod), clockwise=True)

        # print '  get_lens R=%s A=%s N=%s'%(radius, internal_angle_mod, N_turns)

        crow_length = 2 * cl_io + N_turns * cl_wiggle + 2 * inout_space
        path_length = 2 * l_io + N_turns * l_wiggle + 2 * inout_space

        return crow_length, path_length

    # Optimise PATH length with CROW length free
    if target_path_length is not None and target_crow_length is None:
        result = 'path'
        l_targ = target_path_length

    # Optimise CROW length with PATH length free
    elif target_path_length is None and target_crow_length is not None:
        result = 'crow'
        l_targ = target_crow_length
    # Optimise both CROW *AND* PATH
    elif target_path_length is not None and target_crow_length is not None:
        result = 'both'
    # Don't optimise anything, just draw whatever is provided
    else:
        result = 'none'

    # Limits
    r_min = radius
    r_max = 1000.0

    a_min = -89.9 * (2 * np.pi / 360.)
    a_max = 25.0 * (2 * np.pi / 360.)

    if result == 'crow':
        # Function reflecting optimal crow length
        def frootR(radius, internal_angle_mod, N_turns):
            (l_crow, l_path) = get_lens(radius, internal_angle_mod, N_turns)
            return l_crow - target_crow_length

        # Adjust number of bends until the target is in range
        while frootR(r_min, internal_angle_mod, N_turns) > 0 and frootR(r_max, internal_angle_mod, N_turns) > 0:
            N_turns -= 2
        while frootR(r_min, internal_angle_mod, N_turns) < 0 and frootR(r_max, internal_angle_mod, N_turns) < 0:
            N_turns += 2

        # If any bends are required, find the parameter's root
        if N_turns >= 1:
            radius = brentq(frootR, r_min, r_max, args=(internal_angle_mod, N_turns))
        # else:
        #     print('No bends are required')

    elif result == 'path':
        # Function reflecting optimal path length
        def frootA(internal_angle_mod):
            (l_crow, l_path) = get_lens(radius, internal_angle_mod, N_turns)
            return l_path - target_path_length

        # Adjust number of bends until the target is in range
        while frootA(a_min) > 0 and frootA(a_max) > 0:
            N_turns -= 2

        while frootA(a_min) < 0 and frootA(a_max) < 0:
            N_turns += 2

        # If any bends are required, find the parameter's root
        if N_turns >= 1:
            internal_angle_mod = brentq(frootA, a_min, a_max)
        else:
            print('No bends are required')

    elif result == 'both':
        # Function reflecting optimal crow length (optimise radius)
        def frootR(radius, internal_angle_mod, N_turns):
            # print('i_a_m:', internal_angle_mod)
            l_crow, l_path = get_lens(radius, internal_angle_mod, N_turns)
            return l_crow - target_crow_length

        # Wider-scope variables for use inside frootA (below)
        stored_params = [radius, internal_angle_mod, N_turns]

        # Function reflecting optimal path length (optimise internal angle)
        def frootA(internal_angle_mod, radius, N_turns):
            # Adjust number of bends until the target is in range
            while frootR(r_min, internal_angle_mod, N_turns) > 0 and frootR(r_max, internal_angle_mod, N_turns) > 0:
                N_turns -= 2

            while frootR(r_min, internal_angle_mod, N_turns) < 0 and frootR(r_max, internal_angle_mod, N_turns) < 0:
                N_turns += 2

            if frootR(r_min, internal_angle_mod, N_turns) * frootR(r_max, internal_angle_mod, N_turns) > 0:
                raise AttributeError(
                    'Requested geometry (crow length %s, with path length %s) is not possible. Either increase '
                    'target_crow_length or decrease target_path_length.' % (
                        target_crow_length, target_path_length))

            # If any bends are required, find the parameter's root
            if N_turns >= 1:
                radius = brentq(frootR, r_min, r_max, args=(internal_angle_mod, N_turns))
            else:
                print('No bends are required')

            l_crow, l_path = get_lens(radius, internal_angle_mod, N_turns)

            # Store wider-scope variables
            stored_params[0] = radius
            stored_params[1] = internal_angle_mod
            stored_params[2] = N_turns

            return l_path - target_path_length

        # Check that requested geometry is reasonable
        if frootR(r_min, a_max, 1) > 0:
            raise AttributeError(
                'Requested geometry (crow length %s, with radius %s) is not possible.' % (target_crow_length, radius))
        elif target_path_length < target_crow_length:
            raise AttributeError('Target path length (%s) must be greater than target crow length (%s).' % (
                target_path_length, target_crow_length))
        elif target_path_length == target_crow_length:
            N_turns = 0
            l_targ = target_crow_length
        elif frootA(a_min, radius, N_turns) * frootA(a_max, radius, N_turns) > 0:
            raise AttributeError(
                'Requested geometry (crow length %s, with path length %s) is not possible. Either increase '
                'target_crow_length or decrease target_path_length.' % (
                    target_crow_length, target_path_length))
        else:
            # If any bends are required, find the parameter's root
            internal_angle_mod = brentq(frootA, a_min, a_max, args=(radius, N_turns))

            # Restore wider-scope variables
            [radius, internal_angle_mod, N_turns] = stored_params

    ###########
    # DRAWING #
    ###########

    # We move and rotate the shape after drawing it
    wiggle_pts = [np.array((0, 0))]
    if inout_space>0:
        wiggle_pts += [np.array((inout_space, 0))]
    angle = 0

    # If the optimisation returns no bends, shape is a straight line
    if N_turns < 1:
        # Straight line
        wiggle_pts += [np.array((target_crow_length, 0))]

    else:
        # Input bend
        da = (np.pi / 2. + internal_angle_mod)
        temp_pts = EulerBendPoints_Relative(
            start_point=(inout_space, 0),
            radius=radius,
            input_angle=angle,
            angle_amount=da,
            clockwise=False,
            resolution=resolution)
        wiggle_pts += temp_pts[1:]
        angle += da

        # Wiggles
        sign = -1
        for i in range(N_turns):
            clockwise = bool((i + 1) % 2)
            da = sign * (np.pi + 2 * internal_angle_mod)
            temp_pts = EulerBendPoints_Relative(
                start_point=wiggle_pts[-1],
                radius=radius,
                input_angle=angle,
                angle_amount=da,
                clockwise=clockwise,
                resolution=resolution)
            wiggle_pts += temp_pts[1:]
            angle += da
            sign *= -1

        # Output bend
        da = sign * (np.pi + internal_angle_mod)
        temp_pts = EulerBendPoints_Relative(
            start_point=wiggle_pts[-1],
            radius=radius,
            input_angle=angle,
            angle_amount=-angle,
            clockwise=(not clockwise),
            resolution=resolution)
        wiggle_pts += temp_pts[1:]
        angle += da


    last_pt = wiggle_pts[-1]
    if inout_space > 0:
        wiggle_pts += [np.array((last_pt[0] + inout_space, last_pt[1]))]

    # Rotate and mirror shape according to start angle
    if mirrored:
        for pt in wiggle_pts:
            pt[1] *= -1
    origin_pt = np.array(start_point)
    wiggle_pts = [rotate(pt, input_angle, (0, 0)) + origin_pt for pt in wiggle_pts]

    return wiggle_pts


def wgAdd_EulerWiggle(wg,
                      radius=10.0,
                      target_path_length=None,
                      target_crow_length=None,
                      internal_angle_mod=0.0,
                      N_turns=10,
                      mirrored=False,
                      resolution=200):
    """add an EulerWiggle path compensation to a gdshelpers waveguide."""
    wiggle_pts = EulerWiggle_Points(start_point=(0.0, 0.0),
                       input_angle=0.0,
                       radius=radius,
                       target_path_length=target_path_length,
                       target_crow_length=target_crow_length,
                       internal_angle_mod=internal_angle_mod,
                       N_turns=N_turns,
                       mirrored=mirrored,
                       resolution=resolution)
    wiggle_pts = [tuple(pt) for pt in wiggle_pts]
    wg.add_parameterized_path(wiggle_pts)


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from gdshelpers.geometry.chip import Cell
    from gdshelpers.parts.port import Port

    #########################################################
    ##########   Tests bends curves and plot them   #########
    #########################################################

    start_point = (0., 0.)
    radius = 50.0
    input_angle = 0.  # np.pi / 2
    output_angle = np.pi / 3.4
    clockwise = False
    resolution = 200.0

    target_path_length = 510
    target_crow_length = 500
    N_turns = 5
    mirrored = False

    asd = EulerWiggle_Points(start_point=start_point,
                             input_angle=input_angle,
                             radius=radius,
                             target_path_length=target_path_length,
                             target_crow_length=target_crow_length,
                             internal_angle_mod=0.0,
                             N_turns=N_turns,
                             mirrored=mirrored,
                             resolution=resolution)

    # asd = EulerSBend_Points(start_point=start_point,
    #                         offset=-30,
    #                         radius=radius,
    #                         input_angle=input_angle,
    #                         resolution=resolution)

    # asd = EulerBendPoints(start_point=start_point, radius=radius,
    #                       input_angle=input_angle, output_angle=output_angle,
    #                       clockwise=clockwise,
    #                       resolution=resolution)
    #
    # asd_rel = EulerBendPoints_Relative(start_point=start_point, radius=radius,
    #                                    input_angle=input_angle, angle_amount=output_angle - input_angle,
    #                                    clockwise=clockwise,
    #                                    resolution=resolution)
    #
    # end_pt = EulerEndPt(start_point=start_point,
    #                     radius=radius,
    #                     input_angle=input_angle,
    #                     clockwise=clockwise,
    #                     angle_amount=output_angle - input_angle)

    fig, ax = plt.subplots()
    ax.plot([el[0] for el in asd], [el[1] for el in asd], alpha=0.7, lw=3, label='Abs')
    # ax.plot([el[0] for el in asd_rel], [el[1] for el in asd_rel], alpha=0.7, lw=2, label='Rel')
    # print(asd[-1], asd_rel[-1], end_pt)
    ax.legend()
    ax.set_aspect('equal')
    plt.show()

    #########################################################
    ##########       Tests with GDShelpers          #########
    #########################################################

    wg = Waveguide.make_at_port(Port((0, 0), angle=np.pi / 2, width=1.3))
    wg.add_straight_segment(radius)
    wg.add_bend(-np.pi/3, radius, final_width=1.5)

    wgAdd_EulerWiggle(wg,
                      radius=radius,
                      target_path_length=target_path_length,
                      target_crow_length=target_crow_length,
                      internal_angle_mod=0.0,
                      N_turns=N_turns,
                      mirrored=mirrored,
                      resolution=resolution)

    cell = Cell('CELL')
    cell.add_to_layer(1, wg)  # red
    cell.show()