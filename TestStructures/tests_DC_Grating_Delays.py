import numpy as np
import json
from math import pi
import math

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
from gdshelpers.parts.resonator import RingResonator

import os
import sys

p = os.path.abspath('..')
if p not in sys.path:
    sys.path.append(p)

from euler_curves import wgAdd_EulerBend, EulerLength

std_coupler_params = {
    'width': 0.5,
    'full_opening_angle': np.deg2rad(90),
    'grating_period': 0.46,
    'grating_ff': 0.3,
    'n_gratings': 10,
    'ap_max_ff': 0.8,
    'n_ap_gratings': 55,
    'taper_length': 12
}

all_positions = []

wg_width = 0.5
wg_Expwidth = 1.0
l_Exptaper = 40

marker_dims = 20
marker_layer_1 = 3
marker_layer_2 = 4
marker_protection_layer = 15

opt_space = 127
box_dist = 5

euler_to_bend_coeff = 1.8703879865



class Ring_Test:
    def __init__(self, origin, label, gap, ring_r):
        self.origin = origin
        self.label = label

        r_bend = 50
        r_eff = r_bend / euler_to_bend_coeff

        coupler_params = std_coupler_params.copy()

        outports = [Port(self.origin + (opt_space * i, 0), np.pi / 2, wg_width) for i in (0, 1)]
        gratingcouplers = [GratingCoupler.make_traditional_coupler_at_port(outport, **std_coupler_params) for outport in
                           outports]

        port = outports[0].inverted_direction

        wg = Waveguide.make_at_port(port)
        wgAdd_EulerBend(wg, np.pi/2., r_eff, False)

        dist_straight = opt_space - 2 * r_bend
        wg.add_straight_segment(dist_straight/2.)

        ring_res = RingResonator.make_at_port(wg.current_port.inverted_direction, gap=gap, radius=ring_r)

        wg.add_straight_segment(dist_straight/2.)

        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wg.add_straight_segment_until_level_of_port(outports[1])

        label = Text(self.origin + (opt_space/2., 0), 25, self.label, alignment='center-top')

        # self.shapely_object = geometric_union([wg] + [label] + gratingcouplers)
        self.shapely_object = geometric_union(gratingcouplers + [label] + [wg] + [ring_res])

        self.comments = Text(self.origin + (60, 40), 10,
                             'gap={:.2f}\nrad={:.2f}'.format(gap, ring_r),
                             alignment='center-top')


        x_cords = []
        y_cords = []
        box_disty = 2
        for this_poly in self.shapely_object.geoms:
            temp_x_cords, temp_y_cords = this_poly.exterior.coords.xy
            x_cords = x_cords + list(temp_x_cords)
            y_cords = y_cords + list(temp_y_cords)
        x_max = max(x_cords) + box_dist
        x_min = min(x_cords) - box_dist
        y_max = max(y_cords) + box_disty
        y_min = min(y_cords) - box_disty
        box_size = (x_max - x_min, y_max - y_min)

        box = Waveguide(((x_min + x_max) / 2., y_min), np.pi / 2, box_size[0])
        box.add_straight_segment(box_size[1])
        self.box = geometric_union([box])

        self.markers = []
        self.markers_copy = []
        self.markers_protection = []

    def get_shapely_object(self):
        return self.shapely_object

    def get_comments(self):
        return self.comments

    def get_box(self):
        return self.box

    def get_markers(self):
        return self.markers

    def get_markers_copy(self):
        return self.markers_copy

    def get_markers_protection(self):
        return self.markers_protection


class Efficiency_Grating:
    def __init__(self, origin, label, period, ff):
        self.origin = origin
        self.label = label

        r_bend = 50
        r_eff = r_bend / euler_to_bend_coeff

        coupler_params = std_coupler_params.copy()
        coupler_params['grating_period'] = period
        coupler_params['grating_ff'] = ff

        outports = [Port(self.origin + (opt_space * i, 0), np.pi / 2, wg_width) for i in (0, 1)]
        gratingcouplers = [GratingCoupler.make_traditional_coupler_at_port(outport, **std_coupler_params) for outport in
                           outports]

        port = outports[0].inverted_direction

        wg = Waveguide.make_at_port(port)
        wgAdd_EulerBend(wg, np.pi/2., r_eff, False)
        wg.add_straight_segment(opt_space - 2*r_bend)
        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wg.add_straight_segment_until_level_of_port(outports[1])

        label = Text(self.origin + (opt_space/2., 0), 25, self.label, alignment='center-top')

        # self.shapely_object = geometric_union([wg] + [label] + gratingcouplers)
        self.shapely_object = geometric_union(gratingcouplers + [label] + [wg])

        self.comments = Text(self.origin + (60, 40), 10,
                             'period={:.2f}\nff={:.2f}'.format(coupler_params['grating_period'],
                                                               coupler_params['grating_ff']),
                             alignment='center-top')


        x_cords = []
        y_cords = []
        box_disty = 2 #2.5 * marker_dims
        for this_poly in self.shapely_object.geoms:
            temp_x_cords, temp_y_cords = this_poly.exterior.coords.xy
            x_cords = x_cords + list(temp_x_cords)
            y_cords = y_cords + list(temp_y_cords)
        x_max = max(x_cords) + box_dist
        x_min = min(x_cords) - box_dist
        y_max = max(y_cords) + box_disty
        y_min = min(y_cords) - box_disty
        box_size = (x_max - x_min, y_max - y_min)

        box = Waveguide(((x_min + x_max) / 2., y_min), np.pi / 2, box_size[0])
        box.add_straight_segment(box_size[1])
        self.box = geometric_union([box])


        # markers_dist = 50
        # marker_pos = [self.origin + (x_min + 1.2*marker_dims, y_min + 1.2*marker_dims) - self.origin,
        #               self.origin + (x_max - 1.2*marker_dims - markers_dist, y_min + 1.2*marker_dims) - self.origin,
        #               self.origin + (x_min + 1.2*marker_dims + markers_dist, y_max - 1.2*marker_dims) - self.origin]
        # marker_pos1 = [marker_pos[0] + (markers_dist, 0),
        #                marker_pos[1] + (markers_dist, 0),
        #                marker_pos[2] - (markers_dist, 0)]
        #
        # markers = [SquareMarker.make_marker(position, marker_dims)
        #            for position in marker_pos]
        # self.markers = geometric_union(markers)
        #
        # markers1 = [SquareMarker.make_marker(position, marker_dims)
        #             for position in marker_pos1]
        # self.markers_copy = geometric_union(markers1)
        #
        # markers_protection = [SquareMarker.make_marker(position, 2 * marker_dims)
        #                       for position in marker_pos + marker_pos1]
        # self.markers_protection = geometric_union(markers_protection)
        self.markers = []
        self.markers_copy = []
        self.markers_protection = []

    def get_shapely_object(self):
        return self.shapely_object

    def get_comments(self):
        return self.comments

    def get_box(self):
        return self.box

    def get_markers(self):
        return self.markers

    def get_markers_copy(self):
        return self.markers_copy

    def get_markers_protection(self):
        return self.markers_protection


###################################################################################
class DirectionalCouplersTest:
    def __init__(self, origin, label, gap, length, AMZI_DeltaL, r_curve=50):
        self.origin = origin
        self.label = label
        couplerList = []
        outports = []
        wgs = []

        r_eff = r_curve / euler_to_bend_coeff

        port = Port(self.origin, 0, wg_width)

        couplerList.append(DirectionalCoupler.make_at_port(port, length=length, gap=gap, bend_radius=r_curve))
        wg = Waveguide.make_at_port(couplerList[0].right_ports[1])
        wg.add_straight_segment(2 * r_curve)
        wgs.append(wg)

        wg = Waveguide.make_at_port(couplerList[0].right_ports[0])
        wg.add_straight_segment(2 * r_curve)
        couplerList.append(DirectionalCoupler.make_at_port(wg.current_port, length=length, gap=gap, bend_radius=r_curve))

        sep=5
        len_eulers = EulerLength(r_eff, np.pi/2.)
        Delta_L = max(0, (AMZI_DeltaL - 12*len_eulers - 5*sep)/4.)#100


        wg = Waveguide.make_at_port(couplerList[0].right_ports[0])
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wg.add_straight_segment(Delta_L+sep)
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wg.add_straight_segment(Delta_L)
        wgAdd_EulerBend(wg, np.pi/2., r_eff, False)
        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wg.add_straight_segment(Delta_L+sep)
        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wg.add_straight_segment(6*r_curve + sep)
        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wg.add_straight_segment(Delta_L+sep)
        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wgAdd_EulerBend(wg, np.pi / 2., r_eff, False)
        wg.add_straight_segment(Delta_L)
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wg.add_straight_segment(Delta_L+sep)
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wgs.append(wg)
        #
        wg = Waveguide.make_at_port(couplerList[0].left_ports[0])
        wgAdd_EulerBend(wg, -np.pi/2., r_eff, True)
        wgs.append(wg)
        outports.append(wg.port)

        outports.append(outports[0].parallel_offset(-3 * opt_space))

        wg = Waveguide.make_at_port(couplerList[1].right_ports[0])
        wg.add_bezier_to(np.array(outports[1].origin), bend_strength=r_curve / 1.2,
                         final_angle=np.pi / 2.)
        wgs.append(wg)

        wg = Waveguide.make_at_port(couplerList[0].left_ports[1])
        wg.add_bend(angle=-np.pi, radius=r_curve/4)
        wg.add_straight_segment(50, final_width=0.1)
        wgs.append(wg)

        wg = Waveguide.make_at_port(couplerList[1].right_ports[1])
        wg.add_bend(angle=np.pi, radius=r_curve/4)
        wg.add_straight_segment(50, final_width=0.1)
        wgs.append(wg)

        gratingcouplers = [GratingCoupler.make_traditional_coupler_at_port(outport, **std_coupler_params) for outport in
                           outports]

        label = Text(self.origin + (180, 70), 25, self.label, alignment='center-top')

        box = Waveguide(self.origin + (+520, -40), np.pi / 2., 1025)
        box.add_straight_segment(730)

        self.shapely_object = geometric_union(wgs + [label] + couplerList + gratingcouplers)

        self.comments = Text(self.origin + (105, 70), 10, 'gap={:.2f}\nlength={:.2f}'.format(gap, length),
                             alignment='center-top')

        x_cords = []
        y_cords = []
        for this_poly in self.shapely_object.geoms:
            temp_x_cords, temp_y_cords = this_poly.exterior.coords.xy
            x_cords = x_cords + list(temp_x_cords)
            y_cords = y_cords + list(temp_y_cords)
        x_max = max(x_cords) + box_dist
        x_min = min(x_cords) - box_dist
        y_max = max(y_cords) + box_dist
        y_min = min(y_cords) - box_dist
        box_size = (x_max - x_min, y_max - y_min)

        box = Waveguide(((x_min + x_max) / 2., y_min), np.pi / 2, box_size[0])
        box.add_straight_segment(box_size[1])
        self.box = geometric_union([box])

        # markers_dist = 50
        # marker_pos = [self.origin + (-40, -85),
        #               self.origin + (110, -85+markers_dist),
        #               self.origin + (340, -85)]
        # marker_pos1 = [marker_pos[0] + (0, markers_dist),
        #                marker_pos[1] + (0, -markers_dist),
        #                marker_pos[2] + (0, markers_dist)]
        #
        # markers = [SquareMarker.make_marker(position, marker_dims)
        #            for position in marker_pos]
        # self.markers = geometric_union(markers)
        #
        # markers1 = [SquareMarker.make_marker(position, marker_dims)
        #             for position in marker_pos1]
        # self.markers_copy = geometric_union(markers1)
        #
        # markers_protection = [SquareMarker.make_marker(position, 2 * marker_dims)
        #                       for position in marker_pos + marker_pos1]
        # self.markers_protection = geometric_union(markers_protection)
        self.markers = []
        self.markers_copy = []
        self.markers_protection = []

    def get_shapely_object(self):
        return self.shapely_object

    def get_comments(self):
        return self.comments

    def get_box(self):
        return self.box

    def get_couplers(self):
        return self.couplers

    def get_markers(self):
        return self.markers

    def get_markers_copy(self):
        return self.markers_copy

    def get_markers_protection(self):
        return self.markers_protection




##########################################################################
class RectangularSpiral:
    def __init__(self, origin, label, num, add_xlength=0., add_ylength=100., sep=4., r_curve=50.):

        r_eff = r_curve / euler_to_bend_coeff

        delta = 2 * sep
        orizontal_length = add_xlength + 2 * sep
        max_box_size = (1040, 1040)

        self.origin = origin
        self.label = label
        coupler_params = std_coupler_params.copy()

        grating_coupl_pos = self.origin

        couplers = [
            GratingCoupler.make_traditional_coupler(grating_coupl_pos + (opt_space * x, 0), **coupler_params) for x
            in
            (0, 1)]

        ports = [Port(grating_coupl_pos - (opt_space * x, 0), np.pi / 2, wg_width) for x in (0, 1)]

        couplers = [GratingCoupler.make_traditional_coupler_at_port(port, **coupler_params) for port in ports]

        waveguides = [Waveguide.make_at_port(port.inverted_direction) for port in ports]
        waveguides[0].add_straight_segment(sep, final_width=wg_width)
        wgAdd_EulerBend(waveguides[0], -np.pi / 2, r_eff, True)
        waveguides[0].add_straight_segment(sep + opt_space)
        wgAdd_EulerBend(waveguides[0], -np.pi / 2., r_eff, True)

        wgAdd_EulerBend(waveguides[1], -np.pi / 2., r_eff, True)
        wgAdd_EulerBend(waveguides[1], -np.pi / 2., r_eff, True)

        origin_spiral = Port(self.origin + (-opt_space / 2., 320 + add_ylength / 2. + num * sep), 0, wg_width)

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

        waveguides[0].add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        waveguides[0].add_straight_segment_until_y(wg1.y - r_curve - l_Exptaper)
        waveguides[0].add_straight_segment(l_Exptaper, final_width=wg_width)
        wgAdd_EulerBend(waveguides[0], -np.pi / 2., r_eff, True)
        waveguides[0].add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        waveguides[0].add_straight_segment_until_x(wg1.x - l_Exptaper)
        waveguides[0].add_straight_segment(l_Exptaper, final_width=wg_width)

        waveguides[1].add_straight_segment(r_curve / 2., final_width=wg_width)
        wgAdd_EulerBend(waveguides[1], -np.pi / 2., r_eff, True)

        if (num * delta - sep + orizontal_length) > 2 * l_Exptaper:
            wg2.add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
            wg2.add_straight_segment(num * delta - sep + orizontal_length - 2 * l_Exptaper)
            wg2.add_straight_segment(l_Exptaper, final_width=wg_width)
        else:
            wg2.add_straight_segment(num * delta - sep + orizontal_length)
        wgAdd_EulerBend(wg2, -np.pi / 2., r_eff, True)

        waveguides[1].add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        waveguides[1].add_straight_segment_until_x(wg2.x - r_curve - l_Exptaper)
        waveguides[1].add_straight_segment(l_Exptaper, final_width=wg_width)
        wgAdd_EulerBend(waveguides[1], np.pi / 2., r_eff, False)
        waveguides[1].add_straight_segment(l_Exptaper, final_width=wg_Expwidth)
        waveguides[1].add_straight_segment_until_y(wg2.y - l_Exptaper)
        waveguides[1].add_straight_segment(l_Exptaper, final_width=wg_width)

        label = Text(origin_spiral.origin + (0, -r_curve), 25, self.label, alignment='center-top')

        self.shapely_object = geometric_union(waveguides + [wg1, wg2, label] + couplers)

        self.comments = Text(self.origin + (-100, -10), 30,
                             'length={:.2f}'.format(
                                 (wg1.length + wg2.length + waveguides[1].length + waveguides[0].length) * 1e-4),
                             alignment='center-top')

        x_cords = []
        y_cords = []
        for this_poly in self.shapely_object.geoms:
            temp_x_cords, temp_y_cords = this_poly.exterior.coords.xy
            x_cords = x_cords + list(temp_x_cords)
            y_cords = y_cords + list(temp_y_cords)
        x_max = max(x_cords) + box_dist
        x_min = min(x_cords) - box_dist
        y_max = max(y_cords) + box_dist
        y_min = min(y_cords) - box_dist
        box_size = (x_max - x_min, y_max - y_min)

        num_boxes = int(box_size[1] / max_box_size[1]) + 1
        box_x_dim = box_size[0]
        box_y_dim = box_size[1] / num_boxes
        box_list = []
        for i in range(num_boxes):
            box = Waveguide(((x_min + x_max) / 2., y_min + i * box_y_dim), np.pi / 2, box_x_dim)
            box.add_straight_segment(box_y_dim)
            box_list.append(box)
        # self.box = geometric_union(box_list)
        self.box = box_list

        markers_dist = 60
        marker_pos = [self.origin + (-180, +120),
                      self.origin + (x_max - 60, y_min + 175) - self.origin,
                      self.origin + (-opt_space / 2., 800)]
        marker_pos1 = [marker_pos[0] + (0, markers_dist),
                       marker_pos[1] + (0, markers_dist),
                       marker_pos[2] + (0, markers_dist)]

        markers = [SquareMarker.make_marker(position, marker_dims)
                   for position in marker_pos]
        self.markers = geometric_union(markers)

        markers1 = [SquareMarker.make_marker(position, marker_dims)
                    for position in marker_pos1]
        self.markers_copy = geometric_union(markers1)

        markers_protection = [SquareMarker.make_marker(position, 2 * marker_dims)
                              for position in marker_pos + marker_pos1]
        self.markers_protection = geometric_union(markers_protection)

    def get_shapely_object(self):
        return self.shapely_object

    def get_comments(self):
        return self.comments

    def get_box(self):
        return self.box

    def get_couplers(self):
        return self.couplers

    def get_markers(self):
        return self.markers

    def get_markers_copy(self):
        return self.markers_copy

    def get_markers_protection(self):
        return self.markers_protection


#######################################################################################
if __name__ == "__main__":
    wg_layer = 9
    box_layer = 100
    comment_layer = 11
    marker_layer_1 = 3
    marker_layer_2 = 4
    marker_protection_layer = 15

    spacing = np.array((1250, 330))
    spacing2 = np.array((500, 330))

    devices = []

    #### ADD RING TESTS
    devices += [Ring_Test(np.array((420, 250)) * (x, y), id_to_alphanumeric(x, y), gap, ring_r)
                for x, gap in enumerate(np.linspace(0.3, 1, 3))
                for y, ring_r in enumerate(np.linspace(30, 80, 3))]

    #### ADD GRATING COUPLERS TESTS
    # devices += [Efficiency_Grating(np.array((420, 250)) * (x, y) + (5000, 2700), id_to_alphanumeric(x, y), period, ff)
    #             for x, period in enumerate(np.linspace(0.88, 0.92, 3)) for y, ff in
    #             enumerate(np.linspace(0.28, 0.33, 3))]

    #### ADD RINGS TESTS
    # devices += [DeviceCouplerRing(np.array((480, 580)) * (x, y) + (-380, 2700), id_to_alphanumeric(y, x), gap, radius)
    #             for x, gap in enumerate(np.linspace(0.1, 1, 10)) for y, radius in enumerate(np.linspace(50, 200, 6))]

    #### ADD DIRECTIONAL-COUPLERS TESTS

    # def FSR_to_L(FSR_GHz, n_g=1.598562):
    #     # CHECK THAT THE REFRACTIVE INDEX IS CORRECT!!!!
    #     c = 2.998 * 10 ** 8
    #     return c / n_g / FSR_GHz * 10 ** 6 * 10 ** -9
    #
    #
    # dL = FSR_to_L(300)
    #
    # print(dL)
    #
    # spacingx = np.array((1045, 0))
    # spacingy = np.array((0, 750))
    #
    # NumTests_x = 3
    # NumTests_y = 2
    #
    # OriginDCTests = np.array((0, 0))
    #
    # for row, gap in enumerate(np.linspace(0.2, 0.5, NumTests_y)):
    #     for col, length in enumerate(np.linspace(20, 35, NumTests_x)):
    #         devices.append(
    #             DirectionalCouplersTest(OriginDCTests + row * spacingy + col * spacingx, str(row) + ',' + str(col), gap,
    #                                     length, dL, r_curve=35))

    #### ADD SPIRALS
    # init_spirals_pos = np.array((0, 0))
    #
    # scan_list = np.linspace(4, 36, 5)
    # dev_dist_list = (0, 237, 237, 264, 304)
    #
    # temp_pos = init_spirals_pos
    # for x, num in enumerate(scan_list):
    #     this_dev = RectangularSpiral(temp_pos + np.array((dev_dist_list[x], 0)), id_to_alphanumeric(x, 0), 2*int(num/2.),
    #                                  add_xlength=0, add_ylength=4000, sep=5, r_curve=45)
    #     devices += [this_dev]
    #     this_box = this_dev.get_box()
    #     temp_pos = np.array((max((this_box.exterior.coords.xy)[0]), init_spirals_pos[1]))

    ## ADD EVERYTHING
    cell = Cell('Spiralone_LiNb')
    # cell.add(convert_to_gdscad(devices))
    for i, device in enumerate(devices):
        print('adding test',i)
        cell.add_to_layer(wg_layer, device.get_shapely_object())
        device_boxes = device.get_box()
        if isinstance(device_boxes, list):
            for this_box in device.get_box():
                cell.add_to_layer(box_layer, this_box)
        else:
            cell.add_to_layer(box_layer, device_boxes)
        cell.add_to_layer(comment_layer, device.get_comments().get_shapely_object())
        cell.add_to_layer(marker_layer_1, device.get_markers())
        cell.add_to_layer(marker_layer_2, device.get_markers_copy())
        cell.add_to_layer(marker_protection_layer, device.get_markers_protection())

    # print('width: {:.0f}\nheight: {:.0f}'.format(*(cell.bounding_box[1] - cell.bounding_box[0])))

    print('starting device saving')
    cell.save('tests_DC_Grating_Delays.gds')
    print('done saving')