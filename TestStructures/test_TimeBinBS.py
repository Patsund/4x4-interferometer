import numpy as np

from gdshelpers.parts.text import Text
from gdshelpers.parts.waveguide import Waveguide
from gdshelpers.parts.port import Port
from gdshelpers.parts.resonator import RingResonator
# from gdshelpers.parts.splitter import Splitter
from gdshelpers.parts.coupler import GratingCoupler
from gdshelpers.parts.splitter import MMI, DirectionalCoupler
from gdshelpers.parts.marker import SquareMarker
# from gdshelpers.parts.spiral import Spiral
from gdshelpers.geometry.shapely_adapter import geometric_union
from gdshelpers.helpers import id_to_alphanumeric
# import gdsCAD.core
from gdshelpers.geometry.chip import Cell
# from gdshelpers.layout import GridLayout
# from gdshelpers.parts.marker import CrossMarker
from shapely.geometry import Polygon

import os
import sys

p = os.path.abspath('..')
if p not in sys.path:
    sys.path.append(p)

from euler_curves import wgAdd_EulerBend, EulerLength, wgAdd_EulerSBend
from tech.LiNb01 import *


def TimeBin_BS(label, num, add_xlength=0., add_ylength=100., sep=4., r_curve=50., coupler_sep=0.4, coupler_length=30,
               MZ_length=1250, electrodes_sep=1.1, return_xmax=False):
    cell = Cell('TimeBinBS_active' + label)

    x_in = 0
    y_in = 0

    wg_sep = mod_params['wg_sep']

    wg_width_in_mod = wg_Expwidth
    taper_length = l_Exptaper

    x_length_dc = 60

    min_safedist_from_wg = 22
    extra_x_dist_to_spir = min_safedist_from_wg + mod_params['electrode_width']

    ##Generating input and output grating couplers

    coupler_params = std_coupler_params.copy()

    # for j in (0, 1):
    #     incoupler = GratingCoupler.make_traditional_coupler((x_in + j * opt_space, y_in), **coupler_params)
    #     cell.add_to_layer(wg_layer, incoupler)

    # for j in (4, 5):
    #     outcoupler = GratingCoupler.make_traditional_coupler((x_in + j * opt_space, y_in), **coupler_params)
    #     cell.add_to_layer(wg_layer, outcoupler)

    ###Generating waveguides

    mzi_inports = [Port((x_in - j * wg_sep, y_in), -np.pi / 2., wg_width) for j in (0, 1)]
    wg = [Waveguide.make_at_port(inport) for inport in mzi_inports]

    for j in (0, 1):
        ##directional coupler with sinusoidal s-bend

        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length_dc, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length_dc, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length_dc, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length_dc, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        # add taper to multimode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width_in_mod)

        ###reference port for electrodes
        ref_port = wg[j].current_port

        ###straight section
        wg[j].add_straight_segment(MZ_length)

        # add taper to single-mode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width)

        ##directional coupler with sinusoidal s-bend
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length_dc, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length_dc, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length_dc, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length_dc, -np.pi * .5 * np.sin(np.pi * t) * y_length))

    #### Add spiral
    r_eff = r_curve / euler_to_bend_coeff
    delta = 2 * sep
    orizontal_length = add_xlength + 2 * sep

    origin_spiral = Port(mzi_inports[0].origin + (
    r_curve + num * sep + extra_x_dist_to_spir, - add_ylength / 2. - 3 * r_curve - num * sep), 0, wg_width)

    wg1 = Waveguide.make_at_port(origin_spiral)
    wgAdd_EulerBend(wg1, -np.pi / 2., r_eff, True)
    if add_ylength > (4 * l_Exptaper):
        wg1.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg1.add_straight_segment(add_ylength / 2. - 2 * l_Exptaper)
        wg1.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg1.add_straight_segment(add_ylength / 2.)
    wgAdd_EulerBend(wg1, -np.pi / 2., r_eff, True)
    if (add_xlength / 2. + sep) > (2 * l_Exptaper):
        wg1.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg1.add_straight_segment(add_xlength / 2. + sep - 2 * l_Exptaper)
        wg1.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg1.add_straight_segment(add_xlength / 2. + sep)
    wgAdd_EulerBend(wg1, -np.pi / 2., r_eff, True)
    if (sep + add_ylength + 2 * r_curve) > (2 * l_Exptaper):
        wg1.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg1.add_straight_segment(sep + add_ylength + 2 * r_curve - 2 * l_Exptaper)
        wg1.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg1.add_straight_segment(sep + add_ylength + 2 * r_curve)
    wgAdd_EulerBend(wg1, -np.pi / 2., r_eff, True)

    wg2 = Waveguide.make_at_port(origin_spiral.inverted_direction)
    wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)
    if add_ylength > (4 * l_Exptaper):
        wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg2.add_straight_segment(add_ylength / 2. - 2 * l_Exptaper)
        wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg2.add_straight_segment(add_ylength / 2.)
    wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)
    if (add_xlength / 2. + sep) > (2 * l_Exptaper):
        wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg2.add_straight_segment(add_xlength / 2. + sep - 2 * l_Exptaper)
        wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg2.add_straight_segment(add_xlength / 2. + sep)
    wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)
    if (sep + add_ylength + 2 * r_curve) > (2 * l_Exptaper):
        wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg2.add_straight_segment(sep + add_ylength + 2 * r_curve - 2 * l_Exptaper)
        wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg2.add_straight_segment(sep + add_ylength + 2 * r_curve)
    wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)
    if (orizontal_length + sep) > (2 * l_Exptaper):
        wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg2.add_straight_segment(orizontal_length + sep - 2 * l_Exptaper)
        wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg2.add_straight_segment(orizontal_length + sep)
    wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)
    if (sep + delta + add_ylength + 2 * r_curve) > (2 * l_Exptaper):
        wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg2.add_straight_segment(sep + delta + add_ylength + 2 * r_curve - 2 * l_Exptaper)
        wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg2.add_straight_segment(sep + delta + add_ylength + 2 * r_curve)
    wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)

    for i in np.arange(1, num):
        if (orizontal_length + delta * i - sep) > (2 * l_Exptaper):
            wg1.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
            wg1.add_straight_segment(orizontal_length + delta * i - sep - 2 * l_Exptaper)
            wg1.add_straight_segment(l_Exptaper, final_width=wg_width)
        else:
            wg1.add_straight_segment(orizontal_length + delta * i - sep)
        wgAdd_EulerBend(wg1, -np.pi / 2., r_eff, True)
        if (delta * i + sep + add_ylength + 2 * r_curve) > (2 * l_Exptaper):
            wg1.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
            wg1.add_straight_segment(delta * i + sep + add_ylength + 2 * r_curve - 2 * l_Exptaper)
            wg1.add_straight_segment(l_Exptaper, final_width=wg_width)
        else:
            wg1.add_straight_segment(delta * i + sep + add_ylength + 2 * r_curve)
        wgAdd_EulerBend(wg1, -np.pi / 2., r_eff, True)

    for j in np.arange(2, num):
        if (orizontal_length + j * delta - sep) > (2 * l_Exptaper):
            wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
            wg2.add_straight_segment(orizontal_length + j * delta - sep - 2 * l_Exptaper)
            wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
        else:
            wg2.add_straight_segment(orizontal_length + j * delta - sep)
        wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)
        if (delta * (j + 1) - sep + add_ylength + 2 * r_curve) > (2 * l_Exptaper):
            wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
            wg2.add_straight_segment(delta * (j + 1) - sep + add_ylength + 2 * r_curve - 2 * l_Exptaper)
            wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
        else:
            wg2.add_straight_segment(delta * (j + 1) - sep + add_ylength + 2 * r_curve)
        wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)

    #### connect MZI output to delay
    to_port = wg[0].current_port
    wg_to = Waveguide.make_at_port(to_port)
    y_diff_to = abs(to_port.origin[1] - wg1.origin[1])
    if to_port.origin[1] - r_curve < wg1.origin[1]:
        if y_diff_to < 2 * r_curve:
            wg_to.add_straight_segment(2 * r_curve - y_diff_to + 2)
        wgAdd_EulerBend(wg_to, np.pi / 2., r_eff, False)
        wg_to.add_straight_segment(extra_x_dist_to_spir)
        wgAdd_EulerBend(wg_to, np.pi / 2., r_eff, False)
        wgAdd_EulerBend(wg_to, np.pi / 2., r_eff, False)
        wgAdd_EulerBend(wg_to, -np.pi / 2., r_eff, True)
        wg_to.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg_to.add_straight_segment_until_y(wg1.origin[1] - r_curve - l_Exptaper)
        wg_to.add_straight_segment(l_Exptaper, final_width=wg_width)
        wg_to.add_bend(-np.pi / 2, r_curve)
    else:
        if y_diff_to > (2 * l_Exptaper + r_curve):
            wg_to.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
            wg_to.add_straight_segment_until_y(wg1.origin[1] + l_Exptaper + r_curve)
            wg_to.add_straight_segment(l_Exptaper, final_width=wg_width)
        else:
            wg_to.add_straight_segment_until_y(wg1.origin[1] + r_curve)
        wg_to.add_bend(np.pi / 2, r_curve)
    if (wg1.origin[0] - wg_to.current_port.origin[0]) > 2 * l_Exptaper:
        wg_to.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        wg_to.add_straight_segment_until_x(wg1.origin[0] - l_Exptaper)
        wg_to.add_straight_segment(l_Exptaper, final_width=wg_width)
    else:
        wg_to.add_straight_segment_until_x(wg1.origin[0])

    #### connect MZI input to delay
    from_port = mzi_inports[0].inverted_direction
    wg_from = Waveguide.make_at_port(from_port)
    wgAdd_EulerBend(wg_from, -np.pi / 2., r_eff, True)
    wg_from.add_straight_segment_until_x(wg2.origin[0])
    wg_from.add_bend(-np.pi / 2., r_curve)
    wg_from.add_straight_segment_until_y(wg2.origin[1] + r_curve)
    wg_from.add_bend(-np.pi / 2., r_curve)

    io_conns_dist = mod_params['electrode_width'] + min_safedist_from_wg

    #### MZI inconn
    MZIin_port = mzi_inports[1].inverted_direction
    wg_MZIin = Waveguide.make_at_port(MZIin_port)
    wgAdd_EulerBend(wg_MZIin, np.pi / 2., r_eff, False)
    wg_MZIin.add_straight_segment(io_conns_dist + opt_space)
    wgAdd_EulerBend(wg_MZIin, np.pi / 2., r_eff, False)
    wg_MZIin.add_straight_segment(io_conns_dist)

    #### MZI outconn
    MZIout_port = wg[1].current_port
    wg_MZIout = Waveguide.make_at_port(MZIout_port)
    wgAdd_EulerBend(wg_MZIout, -np.pi / 2., r_eff, True)
    wgAdd_EulerBend(wg_MZIout, -np.pi / 2., r_eff, True)
    wgAdd_EulerSBend(wg_MZIout, -mzi_inports[1].origin[0] + wg_MZIout.current_port.origin[0] + io_conns_dist, r_curve)
    wg_MZIout.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
    wg_MZIout.add_straight_segment_until_y(mzi_inports[1].origin[1] - io_conns_dist - l_Exptaper)
    wg_MZIout.add_straight_segment(l_Exptaper, final_width=wg_width)
    wgAdd_EulerBend(wg_MZIout, np.pi / 2., r_eff, False)
    wgAdd_EulerBend(wg_MZIout, np.pi / 2., r_eff, False)

    #### Add In/Out grating couplers
    inoutcouplers = [GratingCoupler.make_traditional_coupler_at_port(this_wg.current_port, **coupler_params)
                     for this_wg in [wg_MZIin, wg_MZIout]]

    cell.add_to_layer(wg_layer, inoutcouplers)

    all_wgs = wg + [wg1, wg2, wg_to, wg_from, wg_MZIin, wg_MZIout]
    for this_wg in all_wgs:
        cell.add_to_layer(wg_layer, this_wg)

    ##MODULATOR ELECTRODES

    electr_width = mod_params['electrode_width']
    sep_econns = mod_params['electrode_sep_y']
    cross_width = mod_params['crossing_width']
    pads_pitch = mod_params['connector_probe_pitch']
    pads_width = mod_params['connector_probe_dims'][0]
    pads_len = mod_params['connector_probe_dims'][1] + 600
    add_dist_wgcross = 70
    add_dist_curvespace = 20
    x_safe_dist = ref_port.origin[0] - min_safedist_from_wg - pads_width - add_dist_wgcross - add_dist_curvespace

    ##left ground electrode
    Inport = Port((ref_port.origin[0] - electrodes_sep / 2.0 - wg_sep / 2.0, ref_port.origin[1]), np.deg2rad(-90),
                  electr_width)
    g_left = Waveguide.make_at_port(Inport)
    g_left.add_straight_segment(MZ_length - cross_width / 2. - 2 * sep_econns)
    g_left._current_port.angle = g_left.current_port.angle - np.pi / 2.0
    g_left._current_port.origin[0] = g_left.current_port.origin[0] + g_left.current_port.width / 2.0
    g_left._current_port.width = cross_width
    g_left.add_straight_segment(add_dist_wgcross)
    g_left.add_straight_segment(wg_sep, electr_width)
    g_left.add_straight_segment_until_x(x_safe_dist)
    g_left.add_straight_segment(2 * pads_pitch)
    g_left._current_port.angle = g_left.current_port.angle + np.pi / 2.0
    g_left._current_port.origin[0] = g_left.current_port.origin[0] + pads_width / 2.0
    g_left._current_port.width = pads_width
    g_left.add_straight_segment(pads_len)
    cell.add_to_layer(electrode_layer, g_left)

    ## signal electrode
    Inport = Port((ref_port.origin[0] + wg_sep / 2.0, ref_port.origin[1]), np.deg2rad(-90),
                  electr_width - electrodes_sep)
    s = Waveguide.make_at_port(Inport)
    s.add_straight_segment(MZ_length - cross_width / 2. - sep_econns)
    s._current_port.angle = s.current_port.angle - np.pi / 2.0
    s._current_port.origin[0] = s.current_port.origin[0] + s.current_port.width / 2.0
    s._current_port.width = cross_width
    s.add_straight_segment(wg_sep + add_dist_wgcross)
    s.add_straight_segment(wg_sep, electr_width)
    s.add_straight_segment_until_x(x_safe_dist)
    s.add_straight_segment(pads_pitch)
    s._current_port.angle = s.current_port.angle + np.pi / 2.0
    s._current_port.origin[0] = s.current_port.origin[0] + pads_width / 2.0
    s._current_port.width = pads_width
    s.add_straight_segment_until_y(g_left.current_port.origin[1] + 50)
    cell.add_to_layer(electrode_layer, s)

    ##right ground electrode
    Inport = Port((ref_port.origin[0] + wg_sep + wg_sep / 2.0 + electrodes_sep / 2.0, ref_port.origin[1]),
                  np.deg2rad(-90), electr_width)
    g_right = Waveguide.make_at_port(Inport)
    g_right.add_straight_segment(MZ_length - cross_width / 2.)
    g_right._current_port.angle = g_right.current_port.angle - np.pi / 2.0
    g_right._current_port.origin[0] = g_right.current_port.origin[0] + g_right.current_port.width / 2.0
    g_right._current_port.width = cross_width
    g_right.add_straight_segment(2 * wg_sep + add_dist_wgcross)
    g_right.add_straight_segment(wg_sep, electr_width)
    g_right.add_straight_segment_until_x(x_safe_dist)
    g_right._current_port.angle = g_right.current_port.angle + np.pi / 2.0
    g_right._current_port.origin[0] = g_right.current_port.origin[0] + pads_width / 2.0
    g_right._current_port.width = pads_width
    g_right.add_straight_segment_until_y(g_left.current_port.origin[1])
    g_right._current_port.angle = g_right.current_port.angle - np.pi / 2.0
    g_right._current_port.origin[0] = g_right.current_port.origin[0] + g_right.current_port.width / 2.0
    g_right._current_port.width = pads_width / 2.
    g_right._current_port.origin[1] = g_right.current_port.origin[1] - g_right._current_port.width / 2.
    g_right.add_straight_segment(2 * pads_pitch + pads_width)
    cell.add_to_layer(electrode_layer, g_right)

    ##WRITE FIELDs waveguide
    added_y_space_formarkers = 160
    x_cords = []
    y_cords = []
    wg_shapely_object = geometric_union(all_wgs + inoutcouplers)
    for this_poly in wg_shapely_object.geoms:
        temp_x_cords, temp_y_cords = this_poly.exterior.coords.xy
        x_cords = x_cords + list(temp_x_cords)
        y_cords = y_cords + list(temp_y_cords)
    x_max = max(x_cords) + box_dist
    x_min = min(x_cords) - box_dist
    y_max = max(y_cords) + box_dist + added_y_space_formarkers
    y_min = min(y_cords) - box_dist
    box_size = (x_max - x_min, y_max - y_min)

    num_boxes = int(box_size[1] / max_box_size[1]) + 1
    box_x_dim = box_size[0]
    box_y_dim = box_size[1] / num_boxes
    box_list = []
    for i in range(num_boxes):
        box = Waveguide(((x_min + x_max) / 2., y_max - i * box_y_dim), -np.pi / 2, box_x_dim)
        box.add_straight_segment(box_y_dim)
        box_list.append(box)
    for this_box in box_list:
        cell.add_to_layer(wg_wf_layer, this_box)
    all_boxes = geometric_union(box_list)
    cell.add_to_layer(wg_reg_layer, all_boxes)

    ##WRITE FIELDs electrodes
    added_x_space_formarkers = 200
    x_cords = []
    y_cords = []
    elec_shapely_object = geometric_union([g_left, s, g_right])
    for this_poly in elec_shapely_object.geoms:
        temp_x_cords, temp_y_cords = this_poly.exterior.coords.xy
        x_cords = x_cords + list(temp_x_cords)
        y_cords = y_cords + list(temp_y_cords)
    x_maxe = min([max(x_cords) + box_dist + added_x_space_formarkers, x_max])
    x_mine = min(x_cords) - box_dist
    y_maxe = y_max  # max(y_cords) + box_dist + 140
    y_mine = min(y_cords) - box_dist
    box_size = (x_maxe - x_mine, y_maxe - y_mine)

    num_boxes = int(box_size[1] / max_box_size[1]) + 1
    box_x_dim = box_size[0]
    box_y_dim = box_size[1] / num_boxes
    box_list = []
    for i in range(num_boxes):
        box = Waveguide(((x_mine + x_maxe) / 2., y_maxe - i * box_y_dim), -np.pi / 2, box_x_dim)
        box.add_straight_segment(box_y_dim)
        box_list.append(box)
    for this_box in box_list:
        cell.add_to_layer(electrode_wf_layer, this_box)
    all_boxes = geometric_union(box_list)
    cell.add_to_layer(electrode_reg_layer, all_boxes)

    ####Local markers
    dist_mtom = 100
    dist_mtobox = 25
    y_dist_save = 420
    # ### first set on layer 3
    positions = [(x_min + dist_mtobox, y_max - y_dist_save - dist_mtobox),
                 (x_min + dist_mtobox, y_max - y_dist_save - dist_mtobox - 2 * dist_mtom),
                 (x_maxe - dist_mtobox, y_max - dist_mtobox)]
    marker = [SquareMarker.make_marker(position, 20) for position in positions]
    cell.add_to_layer(marker_layer_1, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 30) for position in positions]
    cell.add_to_layer(wg_layer, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 40) for position in positions]
    cell.add_to_layer(marker_protection_layer, geometric_union(marker))
    #
    # ### second set on layer 4
    positions = [(x_min + dist_mtobox, y_max - y_dist_save - dist_mtobox - dist_mtom),
                 (x_min + dist_mtobox, y_max - y_dist_save - dist_mtobox - 3 * dist_mtom),
                 (x_maxe - dist_mtobox - dist_mtom, y_max - dist_mtobox)]
    marker = [SquareMarker.make_marker(position, 20) for position in positions]
    cell.add_to_layer(marker_layer_2, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 30) for position in positions]
    cell.add_to_layer(wg_layer, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 40) for position in positions]
    cell.add_to_layer(marker_protection_layer, geometric_union(marker))

    ###Label
    device_label = Text(origin=(x_mine+300, y_max - 50), height=30,
                        text=label, alignment='center-bottom', angle=np.pi)
    cell.add_to_layer(wg_layer, device_label)

    # comments = Text((-100, -10), 30,
    #                 'length={:.2f}'.format(
    #                     (wg1.length + wg2.length + waveguides[1].length + waveguides[0].length) * 1e-4),
    #                 alignment='center-top')
    # cell.add_to_layer(comment_layer, comments)

    ###Device Info
    tot_length = (wg1.length + wg2.length + wg_to.length + wg_from.length) * 1e-4
    info_text = ('delay_length = %.3f cm\nMZ_length = %.1f um\nElectrodes_sep = %.1f nm\nCoupler_length= %.1f um\n') \
                % (tot_length, MZ_length, electrodes_sep, coupler_length)
    device_info = Text(origin=(x_mine+300, y_in - 850), height=20, text=info_text, alignment='center-bottom', angle=np.pi)
    cell.add_to_layer(comment_layer, device_info)

    if return_xmax:
        return cell, -x_mine
    else:
        return cell


#######################################################################################
if __name__ == "__main__":
    devices = []

    #### ADD Modulators
    global_cell = Cell('Modulators_Tests')

    # global_cell.add_to_layer(marker_protection_layer, device.get_markers_protection())

    for j, num in enumerate(np.linspace(10, 30, 2)):
        temp_cell = TimeBin_BS('BS%i' % j, num, add_xlength=0., add_ylength=7000., sep=5., r_curve=50., coupler_sep=0.4,
                               coupler_length=30, MZ_length=1250, electrodes_sep=1.1)
        temp_cell.name = 'TimeBinBS_test_' + str(j)
        global_cell.add_cell(temp_cell, origin=(-2000 + j * 2000, -1500), angle=np.pi)
        # global_cell.add_cell(temp_cell, origin=(-2000 + j * 2000, -1500))

    print('starting device saving')
    global_cell.save('tests_TimeBinBS.gds')
    print('done saving')
