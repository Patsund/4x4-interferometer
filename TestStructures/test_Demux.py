import numpy as np
from copy import deepcopy

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



def Demux_active(coupler_sep, coupler_length, Mod_length, electrodes_sep, label, exp_wg_width=wg_Expwidth, grating_coupler_period = std_coupler_params['grating_period']):
    cell = Cell('Demux_active'+label)

    x_in = 0
    y_in = 0

    MZ_length = Mod_length + 50

    wg_sep = mod_params['wg_sep']
    wg_sep_out = 60
    wg_sep_small = 10

    wg_width_in_mod = wg_Expwidth
    taper_length = l_Exptaper

    ##Generating input and output grating couplers

    coupler_params = std_coupler_params.copy()
    coupler_params['grating_period'] = grating_coupler_period

    for j in range(4):
        incoupler = GratingCoupler.make_traditional_coupler((x_in + j * opt_space, y_in), **coupler_params)
        cell.add_to_layer(wg_layer, incoupler)

    outcouplers = []
    empty_gratcouplers = 0

    for j in range(4):
        outcoupler = GratingCoupler.make_traditional_coupler((x_in + (j + 4 + empty_gratcouplers) * opt_space, y_in), **coupler_params)
        outcouplers.append(outcoupler)
        cell.add_to_layer(wg_layer, outcoupler)

    ###Generating waveguides

    inports = [Port((x_in + j * opt_space, y_in), np.deg2rad(90), std_coupler_params['width']) for j in (0, 1, 2, 3)]
    wg = [Waveguide.make_at_port(inport) for inport in inports]

    x_ref0 = x_in + 1.5 * opt_space
    # y_ref0 = y_in - 3*bend_r

    x_refout = x_ref0 + (4 + empty_gratcouplers)*opt_space

    x_lastout = x_in + (7 + empty_gratcouplers)*opt_space
    x_centre = (x_in + x_lastout)/2.
    y_start_epads = y_in - MZ_length - 930

    ########## Input Connections

    for j in range(4):
        # adding final tapers as suggested by Munster, SP 21/10/21
        wg[j].add_straight_segment(grating_added_taper_len, final_width=wg_width)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment_until_x(x_in-bend_r - 1)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment(grating_added_taper_len + 20)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment_until_x(x_ref0 - (j - 1.5) * wg_sep_out - bend_r)
        wg[j].add_bend(-np.pi / 2.0, bend_r)

    y_ref0 = wg[-1].current_port.origin[1]-1
    for j in range(4):
        wg[j].add_straight_segment_until_y(y_ref0)

    ########## Zero-th MZI
        mzi0_xstart = x_ref0
        mzi0_ystart = y_ref0

    for j in (1, 2):
        wgAdd_EulerSBend(wg[j], offset= (j-1.5) * (wg_sep_out - wg_sep), radius=bend_r)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-1) * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-1) * (wg_sep  - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        #add taper to multimode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width_in_mod)

        ###reference port for electrodes in MZI
        ref_port0 = wg[j].current_port

        ###straight section
        wg[j].add_straight_segment(MZ_length)

        #add taper to single-mode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-1) * (wg_sep  - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-1) * (wg_sep  - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        wgAdd_EulerSBend(wg[j], offset=-(j - 1.5) * (wg_sep_out - wg_sep), radius=bend_r)

    mzi0_xend = x_ref0
    mzi0_yend = wg[j].current_port.origin[1]
    for j in (0, 3):
        wg[j].add_straight_segment(l_Exptaper, final_width=exp_wg_width)
        wg[j].add_straight_segment_until_y(mzi0_yend + l_Exptaper)
        wg[j].add_straight_segment(l_Exptaper, final_width=wg_width)

    ########## Turn
    for j in range(4):
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment(l_Exptaper, final_width=exp_wg_width)
        wg[j].add_straight_segment_until_x(x_refout - 1.5 * wg_sep_out - bend_r + j*(wg_sep_out-wg_sep_small) - l_Exptaper)
        wg[j].add_straight_segment(l_Exptaper, final_width=wg_width)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)

    ########## First MZI (left)
    mzi1_xstart = wg[0].current_port.origin[0]
    mzi1_ystart = wg[0].current_port.origin[1]

    for j in (0, 1):
        wgAdd_EulerSBend(wg[j], offset=(j-0.5) * (wg_sep_out - wg_sep), radius=bend_r)

        #directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        # add taper to multimode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width_in_mod)

        ###reference port for electrodes in MZI
        ref_port1 = wg[j].current_port

        ###straight section
        wg[j].add_straight_segment(MZ_length)

        # add taper to single-mode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        wgAdd_EulerSBend(wg[j], offset=-(j - 0.5) * (wg_sep_out - wg_sep), radius=bend_r)



    ########## Second MZI (left)
    mzi1_xstart = wg[2].current_port.origin[0]
    mzi1_ystart = wg[2].current_port.origin[1]

    for j in (2, 3):
        wgAdd_EulerSBend(wg[j], offset=(j - 2.5) * (wg_sep_out - wg_sep), radius=bend_r)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-2) * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-2) * (wg_sep - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        # add taper to multimode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width_in_mod)

        ###reference port for electrodes in MZI
        ref_port2 = wg[j].current_port

        ###straight section
        wg[j].add_straight_segment(MZ_length)

        # add taper to single-mode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-2) * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - (j-2) * (wg_sep - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (
                                     x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        wgAdd_EulerSBend(wg[j], offset=-(j - 2.5) * (wg_sep_out - wg_sep), radius=bend_r)

    ########## Output Connections

    for j in range(4):
        wg[j].add_straight_segment((3-j)*wg_sep_small+1)
        wg[j].add_bend(-np.pi / 2.0, bend_r)
        wg[j].add_straight_segment_until_x(x_lastout + bend_r + 1)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment(grating_added_taper_len + 20)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment_until_x(x_lastout - j*(opt_space) + bend_r + j * wg_sep_small)
        wg[j].add_bend(np.pi / 2.0, bend_r + j * wg_sep_small)
        wg[j].add_straight_segment(grating_added_taper_len, final_width=std_coupler_params['width'])


        #
    for j in range(4):
        cell.add_to_layer(wg_layer, wg[j])

    ##MODULATORs ELECTRODES

    electr_width = mod_params['electrode_width']
    sep_econns = mod_params['electrode_sep_y']
    cross_width = mod_params['crossing_width']
    pads_pitch = mod_params['connector_probe_pitch']
    pads_width = mod_params['connector_probe_dims'][0]
    pads_width_gnd = mod_params['connector_probe_dims_gnd'][0]
    pads_len = mod_params['connector_probe_dims'][1] - 200
    min_safedist_from_wg = 22

    x_start_epads = x_centre - 3*pads_pitch

    #### ELECTRODES IN ZEROTH MZI
    x_safe_dist = ref_port0.origin[0] - min_safedist_from_wg - pads_width - wg_sep / 2.0 - wg_sep_out

    ##left ground electrode
    Inport = Port((ref_port0.origin[0] - electrodes_sep / 2.0 - wg_sep / 2.0, ref_port0.origin[1]), np.deg2rad(-90), electr_width)
    g_left0 = Waveguide.make_at_port(Inport)
    g_left0.add_straight_segment(Mod_length - cross_width/2. - 2 * sep_econns)
    g_left0._current_port.angle = g_left0.current_port.angle - np.pi / 2.0
    g_left0._current_port.origin[0] = g_left0.current_port.origin[0] + g_left0.current_port.width / 2.0
    g_left0._current_port.width = cross_width
    g_left0.add_straight_segment_until_x(x_safe_dist+wg_sep_out)
    g_left0.add_straight_segment(wg_sep, electr_width)
    g_left0.add_straight_segment_until_x(x_safe_dist)

    g_left0.add_straight_segment(2 * sep_econns)
    g_left0._current_port.angle = g_left0.current_port.angle + np.pi / 2.0
    g_left0._current_port.origin[1] = g_left0.current_port.origin[1] + g_left0.current_port.width / 2.0
    g_left0.add_straight_segment_until_y(y_start_epads - 2 * sep_econns)
    g_left0._current_port.angle = g_left0.current_port.angle + np.pi / 2.0
    g_left0._current_port.origin[0] = g_left0.current_port.origin[0] - g_left0.current_port.width / 2.0
    g_left0._current_port.origin[1] = g_left0.current_port.origin[1] + g_left0.current_port.width / 2.0
    g_left0.add_straight_segment_until_x(x_start_epads)

    g_left0._current_port.angle = g_left0.current_port.angle - np.pi / 2.0
    g_left0._current_port.origin[0] = g_left0.current_port.origin[0] #+ pads_width_gnd / 2.0
    # g_left0._current_port.origin[1] = g_left0.current_port.origin[1] + g_left0.current_port.width / 2.0
    g_left0._current_port.width = pads_width_gnd
    g_left0.add_straight_segment(pads_len)

    cell.add_to_layer(electrode_layer, g_left0)

    ## signal electrode
    Inport = Port((ref_port0.origin[0] + wg_sep / 2.0, ref_port0.origin[1]), np.deg2rad(-90), electr_width - electrodes_sep)
    s0 = Waveguide.make_at_port(Inport)
    s0.add_straight_segment(Mod_length - cross_width/2. - sep_econns)
    s0._current_port.angle = s0.current_port.angle - np.pi / 2.0
    s0._current_port.origin[0] = s0.current_port.origin[0] + s0.current_port.width / 2.0
    s0._current_port.width = cross_width
    s0.add_straight_segment_until_x(x_safe_dist+wg_sep_out)
    s0.add_straight_segment(wg_sep, electr_width)
    s0.add_straight_segment_until_x(x_safe_dist)

    s0.add_straight_segment(sep_econns)
    s0._current_port.angle = s0.current_port.angle + np.pi / 2.0
    s0._current_port.origin[1] = s0.current_port.origin[1] + s0.current_port.width / 2.0
    s0.add_straight_segment_until_y(y_start_epads - sep_econns)
    s0._current_port.angle = s0.current_port.angle + np.pi / 2.0
    s0._current_port.origin[0] = s0.current_port.origin[0] - s0.current_port.width / 2.0
    s0.add_straight_segment_until_x(x_start_epads + pads_pitch)
    #
    s0._current_port.angle = s0.current_port.angle - np.pi / 2.0
    s0._current_port.origin[0] = s0.current_port.origin[0] #+ pads_width / 2.0
    s0._current_port.origin[1] = s0.current_port.origin[1] + s0.current_port.width / 2.0
    s0._current_port.width = pads_width
    s0.add_straight_segment_until_y(g_left0.current_port.origin[1] + 50)

    cell.add_to_layer(electrode_layer, s0)

    # ##right ground electrode
    Inport = Port((ref_port0.origin[0] + wg_sep + wg_sep / 2.0 + electrodes_sep / 2.0, ref_port0.origin[1]), np.deg2rad(-90), electr_width)
    g_right0 = Waveguide.make_at_port(Inport)
    g_right0.add_straight_segment(Mod_length - cross_width/2.)
    g_right0._current_port.angle = g_right0.current_port.angle - np.pi / 2.0
    g_right0._current_port.origin[0] = g_right0.current_port.origin[0] + g_right0.current_port.width / 2.0
    g_right0._current_port.width = cross_width
    g_right0.add_straight_segment_until_x(x_safe_dist+wg_sep_out)
    g_right0.add_straight_segment(wg_sep, electr_width)
    g_right0.add_straight_segment_until_x(x_safe_dist)

    g_right0._current_port.angle = g_right0.current_port.angle + np.pi / 2.0
    g_right0._current_port.origin[1] = g_right0.current_port.origin[1] + g_right0.current_port.width / 2.0
    g_right0.add_straight_segment_until_y(y_start_epads)
    g_right0._current_port.angle = g_right0.current_port.angle + np.pi / 2.0
    g_right0._current_port.origin[0] = g_right0.current_port.origin[0] - g_right0.current_port.width / 2.0
    g_right0.add_straight_segment_until_x(x_start_epads + 2 * pads_pitch)

    g_right0._current_port.angle = g_right0.current_port.angle - np.pi / 2.0
    g_right0._current_port.origin[0] = g_right0.current_port.origin[0] #+ 5. #+ pads_width_gnd / 2.0
    g_right0._current_port.origin[1] = g_right0.current_port.origin[1] + g_right0.current_port.width / 2.0
    g_right0._current_port.width = pads_width_gnd
    g_right0.add_straight_segment_until_y(g_left0.current_port.origin[1])
    g_right0._current_port.angle = g_right0.current_port.angle - np.pi / 2.0
    g_right0._current_port.origin[0] = g_right0.current_port.origin[0] + g_right0.current_port.width / 2.0
    g_right0._current_port.width = pads_width / 2.
    g_right0._current_port.origin[1] = g_right0.current_port.origin[1] - g_right0._current_port.width / 2.
    # g_right0.add_straight_segment(2 * pads_pitch + pads_width)
    g_right0.add_straight_segment(2 * pads_pitch + pads_width_gnd)

    cell.add_to_layer(electrode_layer, g_right0)



    # #### ELECTRODES FOR FIRST MZI
    elec_offset = 90
    x_safe_dist = ref_port1.origin[0] + min_safedist_from_wg + pads_width + wg_sep / 2.0 + 2*wg_sep_out

    ##right ground electrode
    Inport = Port((ref_port1.origin[0] + wg_sep - wg_sep / 2.0 + electrodes_sep / 2.0, ref_port1.origin[1] + MZ_length - elec_offset), np.deg2rad(-90), electr_width)
    g_right1 = Waveguide.make_at_port(Inport)
    g_right1.add_straight_segment(Mod_length - cross_width/2. - 2 * sep_econns)
    g_right1._current_port.angle = g_right1.current_port.angle + np.pi / 2.0
    g_right1._current_port.origin[0] = g_right1.current_port.origin[0] - g_right1.current_port.width / 2.0
    g_right1._current_port.width = cross_width
    g_right1.add_straight_segment_until_x(x_safe_dist-wg_sep_out)
    g_right1.add_straight_segment(wg_sep, electr_width)
    g_right1.add_straight_segment_until_x(x_safe_dist)

    g_right1.add_straight_segment(2*sep_econns)
    g_right1._current_port.angle = g_right1.current_port.angle - np.pi / 2.0
    g_right1._current_port.origin[1] = g_right1.current_port.origin[1] + g_right1.current_port.width / 2.0
    g_right1.add_straight_segment_until_y(y_start_epads - 2*sep_econns)
    g_right1._current_port.angle = g_right1.current_port.angle - np.pi / 2.0
    g_right1._current_port.origin[0] = g_right1.current_port.origin[0] + g_right1.current_port.width / 2.0
    g_right1.add_straight_segment_until_x(x_start_epads + 4 * pads_pitch)

    g_right1._current_port.angle = g_right1.current_port.angle + np.pi / 2.0
    g_right1._current_port.origin[0] = g_right1.current_port.origin[0] #+ 5. #+ pads_width_gnd / 2.0
    g_right1._current_port.origin[1] = g_right1.current_port.origin[1] + g_right1.current_port.width / 2.0
    g_right1._current_port.width = pads_width_gnd
    g_right1.add_straight_segment_until_y(g_left0.current_port.origin[1])
    g_right1._current_port.angle = g_right1.current_port.angle - np.pi / 2.0
    g_right1._current_port.origin[0] = g_right1.current_port.origin[0] + g_right1.current_port.width / 2.0
    g_right1._current_port.width = pads_width / 2.
    g_right1._current_port.origin[1] = g_right1.current_port.origin[1] - g_right1._current_port.width / 2.
    g_right1.add_straight_segment(2 * pads_pitch + pads_width)

    cell.add_to_layer(electrode_layer, g_right1)

    ## signal electrode
    Inport = Port((ref_port1.origin[0] - wg_sep / 2.0, ref_port1.origin[1] + MZ_length - elec_offset), np.deg2rad(-90), electr_width - electrodes_sep)
    s1 = Waveguide.make_at_port(Inport)
    s1.add_straight_segment(Mod_length - cross_width/2. - sep_econns)
    s1._current_port.angle = s1.current_port.angle + np.pi / 2.0
    s1._current_port.origin[0] = s1.current_port.origin[0] - s1.current_port.width / 2.0
    s1._current_port.width = cross_width
    s1.add_straight_segment_until_x(x_safe_dist-wg_sep_out)
    s1.add_straight_segment(wg_sep, electr_width)
    s1.add_straight_segment_until_x(x_safe_dist)

    s1.add_straight_segment(sep_econns)
    s1._current_port.angle = s1.current_port.angle - np.pi / 2.0
    s1._current_port.origin[1] = s1.current_port.origin[1] + s1.current_port.width / 2.0
    s1.add_straight_segment_until_y(y_start_epads - sep_econns)
    s1._current_port.angle = s1.current_port.angle - np.pi / 2.0
    s1._current_port.origin[0] = s1.current_port.origin[0] + s1.current_port.width / 2.0
    s1.add_straight_segment_until_x(x_start_epads + 3 * pads_pitch)

    s1._current_port.angle = s1.current_port.angle + np.pi / 2.0
    s1._current_port.origin[0] = s1.current_port.origin[0] #+ pads_width / 2.0
    s1._current_port.origin[1] = s1.current_port.origin[1] + s1.current_port.width / 2.0
    s1._current_port.width = pads_width
    s1.add_straight_segment_until_y(g_left0.current_port.origin[1] + 50)

    cell.add_to_layer(electrode_layer, s1)

    # ##left ground electrode
    Inport = Port((ref_port1.origin[0] - wg_sep - wg_sep / 2.0 - electrodes_sep / 2.0, ref_port1.origin[1] + MZ_length - elec_offset), np.deg2rad(-90), electr_width)
    g_left1 = Waveguide.make_at_port(Inport)
    g_left1.add_straight_segment(Mod_length - cross_width/2.)
    g_left1._current_port.angle = g_left1.current_port.angle + np.pi / 2.0
    g_left1._current_port.origin[0] = g_left1.current_port.origin[0] - g_left1.current_port.width / 2.0
    g_left1._current_port.width = cross_width
    g_left1.add_straight_segment_until_x(x_safe_dist-wg_sep_out)
    g_left1.add_straight_segment(wg_sep, electr_width)
    g_left1.add_straight_segment_until_x(x_safe_dist)

    g_left1._current_port.angle = g_left1.current_port.angle - np.pi / 2.0
    g_left1._current_port.origin[1] = g_left1.current_port.origin[1] + g_left1.current_port.width / 2.0
    g_left1.add_straight_segment_until_y(y_start_epads)
    g_left1._current_port.angle = g_left1.current_port.angle - np.pi / 2.0
    g_left1._current_port.origin[0] = g_left1.current_port.origin[0] + g_left1.current_port.width / 2.0
    g_left1.add_straight_segment_until_x(x_start_epads + 2 * pads_pitch)

    # g_left1._current_port.angle = g_left1.current_port.angle + np.pi / 2.0
    # g_left1._current_port.origin[0] = g_left1.current_port.origin[0] + pads_width_gnd / 2.0
    # g_left1._current_port.origin[1] = g_left1.current_port.origin[1] + g_left1.current_port.width / 2.0
    # g_left1._current_port.width = pads_width_gnd
    # g_left1.add_straight_segment_until_y(g_left0.current_port.origin[1])

    cell.add_to_layer(electrode_layer, g_left1)


    # #### ELECTRODES FOR SECOND MZI
    # x_safe_dist = x_safe_dist + 3*sep_econns

    ##right ground electrode
    Inport = Port((ref_port2.origin[0] + wg_sep - wg_sep / 2.0 + electrodes_sep / 2.0, ref_port2.origin[1] + elec_offset), np.deg2rad(90), electr_width)
    g_right2 = Waveguide.make_at_port(Inport)
    g_right2.add_straight_segment(Mod_length - cross_width/2. - 2 * sep_econns)
    g_right2._current_port.angle = g_right2.current_port.angle - np.pi / 2.0
    g_right2._current_port.origin[0] = g_right2.current_port.origin[0] - g_right2.current_port.width / 2.0
    g_right2._current_port.width = cross_width
    g_right2.add_straight_segment_until_x(x_safe_dist-wg_sep_out)
    g_right2.add_straight_segment(wg_sep, electr_width)
    g_right2.add_straight_segment_until_x(x_safe_dist)

    g_right2.add_straight_segment(3*sep_econns)
    g_right2._current_port.angle = g_right2.current_port.angle - np.pi / 2.0
    g_right2._current_port.origin[1] = g_right2.current_port.origin[1] + g_right2.current_port.width / 2.0
    g_right2.add_straight_segment_until_y(y_start_epads - 2*sep_econns)
    g_right2._current_port.angle = g_right2.current_port.angle - np.pi / 2.0
    g_right2._current_port.origin[0] = g_right2.current_port.origin[0] + g_right2.current_port.width / 2.0
    g_right2.add_straight_segment_until_x(x_start_epads + 4 * pads_pitch)

    # g_right2._current_port.angle = g_right2.current_port.angle + np.pi / 2.0
    # g_right2._current_port.origin[0] = g_right2.current_port.origin[0] #+ pads_width_gnd / 2.0
    # g_right2._current_port.origin[1] = g_right2.current_port.origin[1] + g_right2.current_port.width / 2.0
    # g_right2._current_port.width = pads_width_gnd
    # g_right2.add_straight_segment_until_y(g_left0.current_port.origin[1])

    cell.add_to_layer(electrode_layer, g_right2)

    ## signal electrode
    Inport = Port((ref_port2.origin[0] - wg_sep / 2.0, ref_port2.origin[1] + elec_offset), np.deg2rad(90), electr_width - electrodes_sep)
    s2 = Waveguide.make_at_port(Inport)
    s2.add_straight_segment(Mod_length - cross_width/2. - sep_econns)
    s2._current_port.angle = s2.current_port.angle - np.pi / 2.0
    s2._current_port.origin[0] = s2.current_port.origin[0] - s2.current_port.width / 2.0
    s2._current_port.width = cross_width
    s2.add_straight_segment_until_x(x_safe_dist-wg_sep_out)
    s2.add_straight_segment(wg_sep, electr_width)
    s2.add_straight_segment_until_x(x_safe_dist)

    s2.add_straight_segment(4*sep_econns)
    s2._current_port.angle = s2.current_port.angle - np.pi / 2.0
    s2._current_port.origin[1] = s2.current_port.origin[1] + s2.current_port.width / 2.0
    s2.add_straight_segment_until_y(y_start_epads - 3*sep_econns)
    s2._current_port.angle = s2.current_port.angle - np.pi / 2.0
    s2._current_port.origin[0] = s2.current_port.origin[0] + s2.current_port.width / 2.0
    s2.add_straight_segment_until_x(x_start_epads + 5 * pads_pitch)
    #
    s2._current_port.angle = s2.current_port.angle + np.pi / 2.0
    s2._current_port.origin[0] = s2.current_port.origin[0] #+ pads_width / 2.0
    s2._current_port.origin[1] = s2.current_port.origin[1] + s2.current_port.width / 2.0
    s2._current_port.width = pads_width
    s2.add_straight_segment_until_y(g_left0.current_port.origin[1] + 50)

    cell.add_to_layer(electrode_layer, s2)

    # ##left ground electrode
    Inport = Port((ref_port2.origin[0] - wg_sep - wg_sep / 2.0 - electrodes_sep / 2.0, ref_port2.origin[1] + elec_offset), np.deg2rad(90), electr_width)
    g_left2 = Waveguide.make_at_port(Inport)
    g_left2.add_straight_segment(Mod_length - cross_width/2.)
    g_left2._current_port.angle = g_left2.current_port.angle - np.pi / 2.0
    g_left2._current_port.origin[0] = g_left2.current_port.origin[0] - g_left2.current_port.width / 2.0
    g_left2._current_port.width = cross_width
    g_left2.add_straight_segment_until_x(x_safe_dist-wg_sep_out)
    g_left2.add_straight_segment(wg_sep, electr_width)
    g_left2.add_straight_segment_until_x(x_safe_dist)

    g_left2.add_straight_segment(5*sep_econns)
    g_left2._current_port.angle = g_left2.current_port.angle - np.pi / 2.0
    g_left2._current_port.origin[1] = g_left2.current_port.origin[1] + g_left2.current_port.width / 2.0
    g_left2.add_straight_segment_until_y(y_start_epads - 4*sep_econns)
    g_left2._current_port.angle = g_left2.current_port.angle - np.pi / 2.0
    g_left2._current_port.origin[0] = g_left2.current_port.origin[0] + g_left2.current_port.width / 2.0
    # g_left2._current_port.origin[1] = g_left2.current_port.origin[1] + g_left2.current_port.width / 2.0
    g_left2.add_straight_segment_until_x(x_start_epads + 6 * pads_pitch)
    g_left2._current_port.origin[1] = g_left2.current_port.origin[1] - g_left2.current_port.width / 2.0

    g_left2._current_port.angle = g_left2.current_port.angle + np.pi / 2.0
    g_left2._current_port.origin[0] = g_left2.current_port.origin[0] #+ 5. #+ pads_width_gnd / 2.0
    g_left2._current_port.origin[1] = g_left2.current_port.origin[1] + g_left2.current_port.width
    g_left2._current_port.width = pads_width_gnd
    g_left2.add_straight_segment_until_y(g_left0.current_port.origin[1])
    g_left2._current_port.angle = g_left2.current_port.angle - np.pi / 2.0
    g_left2._current_port.origin[0] = g_left2.current_port.origin[0] + g_left2.current_port.width / 2.0
    g_left2._current_port.width = pads_width / 2.
    g_left2._current_port.origin[1] = g_left2.current_port.origin[1] - g_left2._current_port.width / 2.
    g_left2.add_straight_segment(2 * pads_pitch + pads_width)

    cell.add_to_layer(electrode_layer, g_left2)


    #
    #
    #
    # ###WRITE FIELDs waveguide
    #
    # outer_corners = [(x_in - 80, y_in + 160), (x_in + 5 * 127 + 60, y_in + 160),
    #                  (x_in + 5 * 127 + 60, y_in + 160 - 1040), (x_in - 80, y_in + 160 - 1040)]
    # polygon1 = Polygon(outer_corners)
    # cell.add_to_layer(wg_wf_layer, polygon1)
    # outer_corners = [(x_in - 80, y_in + 160 - 1040), (x_in + 5 * 127 + 60, y_in + 160 - 1040),
    #                  (x_in + 5 * 127 + 60, y_in + 160 - 1040 - (MZ_length + 2*taper_length + 660 - 1040)),
    #                  (x_in - 80, y_in + 160 - 1040 - ((MZ_length + 2*taper_length + 660 - 1040)))]
    # polygon2 = Polygon(outer_corners)
    # cell.add_to_layer(wg_wf_layer, polygon2)
    # polygon = geometric_union([polygon1, polygon2])
    # cell.add_to_layer(wg_reg_layer, polygon)
    #
    # ###WRITE FIELDs electrodes
    #
    # outer_corners = [(x_in - 210, y_in - 100), (x_in + 5 * 127 + 140, y_in - 100),
    #                  (x_in + 5 * 127 + 140, y_in - 100 - 1040), (x_in - 210, y_in - 100 - 1040)]
    # polygon1 = Polygon(outer_corners)
    # cell.add_to_layer(electrode_wf_layer, polygon1)
    # outer_corners = [(x_in - 210, y_in - 100 - 1040), (x_in + 5 * 127 + 140, y_in - 100 - 1040),
    #                  (x_in + 5 * 127 + 140, y_in - 100 - 1040 - (MZ_length + 2*taper_length + 762 - 1040)),
    #                  (x_in - 210, y_in - 100 - 1040 - (MZ_length + 2*taper_length + 762 - 1040))]
    # polygon2 = Polygon(outer_corners)
    # cell.add_to_layer(electrode_wf_layer, polygon2)
    # polygon = geometric_union([polygon1, polygon2])
    # cell.add_to_layer(electrode_reg_layer, polygon)
    #
    # ####Local markers
    #
    # ### first set on layer 3
    # positions = [(x_in, y_in - 320), (x_in + 5 * 127 - 60, y_in - 320), (x_in + 5 * 127 - 60, y_in - 320 - 450)]
    # marker = [SquareMarker.make_marker(position, 20) for position in positions]
    # cell.add_to_layer(3, geometric_union(marker))
    # marker = [SquareMarker.make_marker(position, 30) for position in positions]
    # cell.add_to_layer(9, geometric_union(marker))
    # marker = [SquareMarker.make_marker(position, 40) for position in positions]
    # cell.add_to_layer(15, geometric_union(marker))
    #
    # ### second set on layer 4
    # positions = [(x_in, y_in - 320 - 150), (x_in + 5 * 127 - 60, y_in - 320 - 150), (x_in + 5 * 127-60, y_in - 320 - 300)]
    # marker = [SquareMarker.make_marker(position, 20) for position in positions]
    # cell.add_to_layer(marker_layer_1, geometric_union(marker))
    # marker = [SquareMarker.make_marker(position, 30) for position in positions]
    # cell.add_to_layer(wg_layer, geometric_union(marker))
    # marker = [SquareMarker.make_marker(position, 40) for position in positions]
    # cell.add_to_layer(marker_protection_layer, geometric_union(marker))

    ###Label
    device_label = Text(origin=(x_in, y_in - 700), height=30,
                        text=label, alignment='center-bottom', angle=np.pi)
    cell.add_to_layer(wg_layer, device_label)

    ###Device Info
    info_text = ('Mod_length = %.1f um\nElectrodes_sep = %.1f nm\nCoupler_length= %.1f um\n') \
                % (Mod_length, electrodes_sep, coupler_length)
    device_info = Text(origin=(x_in, y_in - 850), height=20, text=info_text, alignment='center-bottom' , angle=np.pi)
    cell.add_to_layer(comment_layer, device_info)

    return cell




#######################################################################################
if __name__ == "__main__":
    devices = []

    #### ADD Modulators
    global_cell = Cell('Demux_Tests')


    # global_cell.add_to_layer(marker_protection_layer, device.get_markers_protection())

    coupler_sep = 0.5
    coupler_length = 30
    electrodes_sep = 1.1

    for j, MZ_length in enumerate(np.linspace(1250, 1500, 2)):
        temp_cell = Demux_active(coupler_sep, coupler_length, MZ_length, electrodes_sep, 'D%i' % j)
        temp_cell.name = 'Demux_test_' + str(j)
        # global_cell.add_cell(temp_cell, origin=(-2000 + j * 1000, -1500), angle=np.pi)
        global_cell.add_cell(temp_cell, origin=(-2000 + j * 3000, -1500))


    print('starting device saving')
    global_cell.save('tests_Demux.gds')
    print('done saving')