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

from euler_curves import wgAdd_EulerBend, EulerLength
from tech.LiNb01 import *



def MZI_active_with_phase(coupler_sep, coupler_length, MZ_length, electrodes_sep, label):
    cell = Cell('MZI_active_withphase'+label)

    x_in = 0
    y_in = 0

    wg_sep = mod_params['wg_sep']

    wg_width_in_mod = wg_Expwidth
    taper_length = l_Exptaper

    ##Generating input and output grating couplers

    coupler_params = std_coupler_params.copy()

    for j in (0, 1):
        incoupler = GratingCoupler.make_traditional_coupler((x_in + j * opt_space, y_in), **coupler_params)
        cell.add_to_layer(wg_layer, incoupler)

    outcouplers = []
    for j in (3, 4):
        outcoupler = GratingCoupler.make_traditional_coupler((x_in + j * opt_space, y_in), **coupler_params)
        outcouplers.append(outcoupler)
        cell.add_to_layer(wg_layer, outcoupler)

    ###Generating waveguides

    inports = [Port((x_in + j * opt_space, y_in), np.deg2rad(90), std_coupler_params['width']) for j in (0, 1)]
    wg = [Waveguide.make_at_port(inport) for inport in inports]

    for j in (0, 1):
        # adding final tapers as suggested by Munster, SP 21/10/21
        wg[j].add_straight_segment(grating_added_taper_len, final_width=wg_width)
        # wg[j].add_straight_segment(bend_r-grating_added_taper_len)

        wg[j].add_bend(-np.pi / 2.0, bend_r + (1 - j) * wg_sep)
        wg[j].add_straight_segment((1 - j) * (opt_space - wg_sep))
        wg[j].add_bend(-np.pi / 2.0, bend_r + (1 - j) * wg_sep)

        wg[j].add_straight_segment(bend_r)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep  - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))

        #add taper to multimode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width_in_mod)

        ###reference port for electrodes in MZI
        ref_port = wg[j].current_port

        ###reference port for electrodes in PHASE
        ref_port_ph = wg[j].current_port

        ###straight section
        wg[j].add_straight_segment(MZ_length)

        #add taper to single-mode wg
        wg[j].add_straight_segment(taper_length, final_width=wg_width)

        ##directional coupler with sinusoidal s-bend
        x_length = 60
        y_length = wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep  - wg_width - coupler_sep)
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))
        wg[j].add_straight_segment(coupler_length)
        y_length = -(wg_sep / 2.0 - (coupler_sep + wg_width) / 2.0 - j * (wg_sep  - wg_width - coupler_sep))
        wg[j].add_parameterized_path(path=lambda t: (t * x_length, .5 * (np.cos(np.pi * t) - 1) * y_length),
                                     path_derivative=lambda t: (x_length, -np.pi * .5 * np.sin(np.pi * t) * y_length))



        x_dist_inout=2*mod_params['electrode_width'] + 15
        wg[j].add_straight_segment(10)
        wg[j].add_bend(np.pi/2., bend_r + j * wg_sep)
        wg[j].add_straight_segment(x_dist_inout+wg_sep)
        wg[j].add_bend(np.pi, bend_r + j * wg_sep)
        wg[j].add_bend(-np.pi / 2., bend_r + (1-j) * wg_sep)
        wg[j].add_straight_segment_until_y(inports[0].y - 3*bend_r)
        wg[j].add_bend(-np.pi / 2., bend_r + (1 - j) * wg_sep)
        wg[j].add_straight_segment_until_x(outcouplers[-1].origin[0]+bend_r)
        wg[j].add_bend(np.pi / 2., bend_r + j * wg_sep)

        wg[j].add_straight_segment_until_y(inports[0].y)
        wg[j].add_straight_segment(grating_added_taper_len)
        wg[j].add_bend(np.pi / 2., bend_r + j* wg_sep)
        wg[j].add_straight_segment_until_x(outcouplers[1-j].origin[0]+(bend_r + j * wg_sep))
        wg[j].add_bend(np.pi / 2., bend_r + j*wg_sep)
        # # adding final tapers as suggested by Munster, SP 21/10/21
        wg[j].add_straight_segment(grating_added_taper_len, final_width=std_coupler_params['width'])

    for j in (0, 1):
        cell.add_to_layer(wg_layer, wg[j])

    ##MODULATOR ELECTRODES


    #### ELECTRODES IN MZI

    electr_width = mod_params['electrode_width']
    sep_econns = mod_params['electrode_sep_y']
    cross_width = mod_params['crossing_width']
    pads_pitch = mod_params['connector_probe_pitch']
    pads_width = mod_params['connector_probe_dims'][0]
    pads_width_gnd = mod_params['connector_probe_dims_gnd'][0]
    pads_len = mod_params['connector_probe_dims'][1]
    min_safedist_from_wg = 22
    x_safe_dist = ref_port.origin[0] - min_safedist_from_wg - pads_width

    ##left ground electrode
    Inport = Port((ref_port.origin[0] - electrodes_sep / 2.0 - wg_sep / 2.0, ref_port.origin[1]), np.deg2rad(-90), electr_width)
    g_left = Waveguide.make_at_port(Inport)
    g_left.add_straight_segment(MZ_length - cross_width/2. - 2 * sep_econns)
    g_left._current_port.angle = g_left.current_port.angle - np.pi / 2.0
    g_left._current_port.origin[0] = g_left.current_port.origin[0] + g_left.current_port.width / 2.0
    g_left.add_straight_segment_until_x(x_safe_dist)
    g_left.add_straight_segment(2 * pads_pitch + (pads_width-pads_width_gnd/2.))
    g_left._current_port.angle = g_left.current_port.angle + np.pi / 2.0
    g_left._current_port.origin[0] = g_left.current_port.origin[0] + pads_width_gnd / 2.0
    g_left._current_port.width = pads_width_gnd

    g_left.add_straight_segment(pads_len)
    cell.add_to_layer(electrode_layer, g_left)

    ## signal electrode
    Inport = Port((ref_port.origin[0] + wg_sep / 2.0, ref_port.origin[1]), np.deg2rad(-90), electr_width - electrodes_sep)
    s = Waveguide.make_at_port(Inport)
    s.add_straight_segment(MZ_length - cross_width/2. - sep_econns)
    s._current_port.angle = s.current_port.angle - np.pi / 2.0
    s._current_port.origin[0] = s.current_port.origin[0] + s.current_port.width / 2.0
    s._current_port.width = cross_width
    s.add_straight_segment(wg_sep+5)
    s.add_straight_segment(wg_sep, electr_width)
    s.add_straight_segment_until_x(x_safe_dist)
    s.add_straight_segment(pads_pitch+ (pads_width-pads_width_gnd/2.))
    s._current_port.angle = s.current_port.angle + np.pi / 2.0
    s._current_port.origin[0] = s.current_port.origin[0] + pads_width / 2.0
    s._current_port.width = pads_width

    s.add_straight_segment_until_y(g_left.current_port.origin[1] + 50)
    cell.add_to_layer(electrode_layer, s)

    ##right ground electrode
    Inport = Port((ref_port.origin[0] + wg_sep + wg_sep / 2.0 + electrodes_sep / 2.0, ref_port.origin[1]), np.deg2rad(-90), electr_width)
    g_right = Waveguide.make_at_port(Inport)
    g_right.add_straight_segment(MZ_length - cross_width/2.)
    g_right._current_port.angle = g_right.current_port.angle - np.pi / 2.0
    g_right._current_port.origin[0] = g_right.current_port.origin[0] + g_right.current_port.width / 2.0
    g_right._current_port.width = cross_width
    g_right.add_straight_segment(2*wg_sep+5)
    g_right.add_straight_segment(wg_sep, electr_width)
    g_right.add_straight_segment_until_x(x_safe_dist)
    g_right.add_straight_segment(pads_width_gnd/2.)
    g_right._current_port.angle = g_right.current_port.angle + np.pi / 2.0
    g_right._current_port.origin[0] = g_right.current_port.origin[0] + pads_width_gnd / 2.0
    g_right._current_port.width = pads_width_gnd
    g_right.add_straight_segment_until_y(g_left.current_port.origin[1])
    g_right._current_port.angle = g_right.current_port.angle - np.pi / 2.0
    g_right._current_port.origin[0] = g_right.current_port.origin[0] + g_right.current_port.width / 2.0
    g_right._current_port.width = pads_width / 2.
    g_right._current_port.origin[1] = g_right.current_port.origin[1] - g_right._current_port.width / 2.
    g_right.add_straight_segment(2 * pads_pitch + pads_width)

    cell.add_to_layer(electrode_layer, g_right)



    #### ELECTRODES FOR PHASE
    ref_port_ph = deepcopy(ref_port)
    ref_port_ph.origin[0] = ref_port.origin[0] + 2*wg_sep + x_dist_inout
    x_safe_dist_ph = ref_port_ph.origin[0] + min_safedist_from_wg + pads_width
    x_startpos_pads_ph = s._current_port.origin[0] + 4*pads_pitch + (pads_width_gnd - pads_width)/2. + (pads_width-pads_width_gnd/2.)
    y_startpos_pads_ph = s._current_port.origin[1] +280

    ##right ground electrode
    Inport = Port((ref_port_ph.origin[0] + electrodes_sep / 2.0 + wg_sep / 2.0, ref_port.origin[1]), np.deg2rad(-90), electr_width)
    g_right_ph = Waveguide.make_at_port(Inport)
    g_right_ph.add_straight_segment(MZ_length - cross_width/2. - 2 * sep_econns)
    g_right_ph._current_port.angle = g_right_ph.current_port.angle + np.pi / 2.0
    g_right_ph._current_port.origin[0] = g_right_ph.current_port.origin[0] - g_right_ph.current_port.width / 2.0
    g_right_ph.add_straight_segment_until_x(x_startpos_pads_ph)
    g_right_ph._current_port.angle = g_right_ph.current_port.angle - np.pi / 2.0
    g_right_ph._current_port.origin[0] = g_right_ph.current_port.origin[0] - g_right_ph.current_port.width / 2.0
    g_right_ph.add_straight_segment_until_y(y_startpos_pads_ph-2*sep_econns+ g_right_ph.current_port.width)
    g_right_ph._current_port.origin[0] = g_right_ph._current_port.origin[0] - (pads_width_gnd- g_right_ph.current_port.width) / 2.0
    g_right_ph._current_port.width = pads_width_gnd
    g_right_ph.add_straight_segment_until_y(g_left._current_port.origin[1])

    cell.add_to_layer(electrode_layer, g_right_ph)

    ## signal electrode
    Inport = Port((ref_port_ph.origin[0] - wg_sep / 2.0, ref_port_ph.origin[1]), np.deg2rad(-90), electr_width - electrodes_sep)
    s_ph = Waveguide.make_at_port(Inport)
    s_ph.add_straight_segment(MZ_length - cross_width/2. - sep_econns)
    s_ph._current_port.angle = s.current_port.angle + np.pi / 2.0
    s_ph._current_port.origin[0] = s_ph.current_port.origin[0] - s_ph.current_port.width / 2.0
    s_ph._current_port.width = cross_width
    s_ph.add_straight_segment(wg_sep+5)
    s_ph.add_straight_segment(wg_sep, electr_width)
    s_ph.add_straight_segment_until_x(x_startpos_pads_ph - sep_econns)
    s_ph._current_port.angle = s_ph.current_port.angle - np.pi / 2.0
    s_ph._current_port.origin[0] = s_ph.current_port.origin[0] - s_ph.current_port.width / 2.0
    s_ph.add_straight_segment_until_y(y_startpos_pads_ph-sep_econns)
    s_ph._current_port.angle = s_ph.current_port.angle - np.pi / 2.0
    s_ph._current_port.origin[1] = s_ph.current_port.origin[1] + s_ph.current_port.width / 2.0

    s_ph.add_straight_segment_until_x(s.current_port.origin[0] + 3*pads_pitch-pads_width / 2.0)
    s_ph._current_port.angle = s_ph.current_port.angle + np.pi / 2.0
    s_ph._current_port.origin[0] = s_ph.current_port.origin[0] + pads_width / 2.0
    s_ph._current_port.width = pads_width
    s_ph.add_straight_segment_until_y(s.current_port.origin[1])
    cell.add_to_layer(electrode_layer, s_ph)

    ##left ground electrode
    Inport = Port((ref_port_ph.origin[0] - wg_sep - wg_sep / 2.0 - electrodes_sep / 2.0, ref_port_ph.origin[1]), np.deg2rad(-90), electr_width)
    g_left_ph = Waveguide.make_at_port(Inport)
    g_left_ph.add_straight_segment(MZ_length - cross_width/2.)
    g_left_ph._current_port.angle = g_left_ph.current_port.angle + np.pi / 2.0
    g_left_ph._current_port.origin[0] = g_left_ph.current_port.origin[0] - g_left_ph.current_port.width / 2.0
    g_left_ph._current_port.width = cross_width
    g_left_ph.add_straight_segment(2*wg_sep+5)
    g_left_ph.add_straight_segment(wg_sep, electr_width)
    g_left_ph.add_straight_segment_until_x(x_startpos_pads_ph - 2*sep_econns)
    g_left_ph._current_port.angle = g_left_ph.current_port.angle - np.pi / 2.0
    g_left_ph._current_port.origin[0] = g_left_ph.current_port.origin[0] - g_left_ph.current_port.width / 2.0
    g_left_ph.add_straight_segment_until_y(y_startpos_pads_ph)

    g_left_ph._current_port.angle = g_left_ph.current_port.angle - np.pi / 2.0
    g_left_ph._current_port.origin[1] = g_left_ph.current_port.origin[1] + g_left_ph.current_port.width / 2.0
    g_left_ph.add_straight_segment_until_x(s.current_port.origin[0] + 2*pads_pitch - pads_width / 2.0)
    g_left_ph._current_port.angle = g_left_ph.current_port.angle + np.pi / 2.0
    g_left_ph._current_port.origin[0] = g_left_ph.current_port.origin[0] + pads_width_gnd / 2.0
    g_left_ph._current_port.width = pads_width_gnd
    g_left_ph.add_straight_segment_until_y(g_right_ph.current_port.origin[1])

    g_left_ph._current_port.angle = g_left_ph.current_port.angle + np.pi / 2.0
    g_left_ph._current_port.origin[0] = g_left_ph.current_port.origin[0] - g_left_ph.current_port.width / 2.0
    g_left_ph._current_port.width = pads_width / 2.
    g_left_ph._current_port.origin[1] = g_left_ph.current_port.origin[1] - g_left_ph._current_port.width / 2.
    g_left_ph.add_straight_segment(2 * pads_pitch + pads_width)

    cell.add_to_layer(electrode_layer, g_left_ph)



    ###WRITE FIELDs waveguide

    outer_corners = [(x_in - 80, y_in + 160), (x_in + 5 * 127 + 60, y_in + 160),
                     (x_in + 5 * 127 + 60, y_in + 160 - 1040), (x_in - 80, y_in + 160 - 1040)]
    polygon1 = Polygon(outer_corners)
    cell.add_to_layer(wg_wf_layer, polygon1)
    outer_corners = [(x_in - 80, y_in + 160 - 1040), (x_in + 5 * 127 + 60, y_in + 160 - 1040),
                     (x_in + 5 * 127 + 60, y_in + 160 - 1040 - (MZ_length + 2*taper_length + 660 - 1040)),
                     (x_in - 80, y_in + 160 - 1040 - ((MZ_length + 2*taper_length + 660 - 1040)))]
    polygon2 = Polygon(outer_corners)
    cell.add_to_layer(wg_wf_layer, polygon2)
    polygon = geometric_union([polygon1, polygon2])
    cell.add_to_layer(wg_reg_layer, polygon)

    ###WRITE FIELDs electrodes

    outer_corners = [(x_in - 210, y_in - 100), (x_in + 5 * 127 + 140, y_in - 100),
                     (x_in + 5 * 127 + 140, y_in - 100 - 1040), (x_in - 210, y_in - 100 - 1040)]
    polygon1 = Polygon(outer_corners)
    cell.add_to_layer(electrode_wf_layer, polygon1)
    outer_corners = [(x_in - 210, y_in - 100 - 1040), (x_in + 5 * 127 + 140, y_in - 100 - 1040),
                     (x_in + 5 * 127 + 140, y_in - 100 - 1040 - (MZ_length + 2*taper_length + 762 - 1040)),
                     (x_in - 210, y_in - 100 - 1040 - (MZ_length + 2*taper_length + 762 - 1040))]
    polygon2 = Polygon(outer_corners)
    cell.add_to_layer(electrode_wf_layer, polygon2)
    polygon = geometric_union([polygon1, polygon2])
    cell.add_to_layer(electrode_reg_layer, polygon)

    ####Local markers

    ### first set on layer 3
    positions = [(x_in, y_in - 320), (x_in + 5 * 127 - 60, y_in - 320), (x_in + 5 * 127 - 60, y_in - 320 - 450)]
    marker = [SquareMarker.make_marker(position, 20) for position in positions]
    cell.add_to_layer(3, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 30) for position in positions]
    cell.add_to_layer(9, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 40) for position in positions]
    cell.add_to_layer(15, geometric_union(marker))

    ### second set on layer 4
    positions = [(x_in, y_in - 320 - 150), (x_in + 5 * 127 - 60, y_in - 320 - 150), (x_in + 5 * 127-60, y_in - 320 - 300)]
    marker = [SquareMarker.make_marker(position, 20) for position in positions]
    cell.add_to_layer(marker_layer_1, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 30) for position in positions]
    cell.add_to_layer(wg_layer, geometric_union(marker))
    marker = [SquareMarker.make_marker(position, 40) for position in positions]
    cell.add_to_layer(marker_protection_layer, geometric_union(marker))

    ###Label
    device_label = Text(origin=(x_in, y_in - 700), height=30,
                        text=label, alignment='center-bottom', angle=np.pi)
    cell.add_to_layer(wg_layer, device_label)

    ###Device Info
    info_text = ('MZ_length = %.1f um\nElectrodes_sep = %.1f nm\nCoupler_length= %.1f um\n') \
                % (MZ_length, electrodes_sep, coupler_length)
    device_info = Text(origin=(x_in, y_in - 850), height=20, text=info_text, alignment='center-bottom' , angle=np.pi)
    cell.add_to_layer(comment_layer, device_info)

    return cell




#######################################################################################
if __name__ == "__main__":
    devices = []

    #### ADD Modulators
    global_cell = Cell('Modulators_Tests')


    # global_cell.add_to_layer(marker_protection_layer, device.get_markers_protection())

    coupler_sep = 0.4
    coupler_length = 22
    electrodes_sep = 1.1

    for j, MZ_length in enumerate(np.linspace(500, 1250, 4)):
        temp_cell = MZI_active_with_phase(coupler_sep, coupler_length, MZ_length, electrodes_sep, 'D%i' % j)
        temp_cell.name = 'MZI_test_' + str(j)
        # global_cell.add_cell(temp_cell, origin=(-2000 + j * 1000, -1500), angle=np.pi)
        global_cell.add_cell(temp_cell, origin=(-2000 + j * 1000, -1500))


    print('starting device saving')
    global_cell.save('tests_ModulatorsWithPhase.gds')
    print('done saving')