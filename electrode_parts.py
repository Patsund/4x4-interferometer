import numpy as np
from utils import *
from gdshelpers.parts.waveguide import Waveguide
from gdshelpers.parts.port import Port
from gdshelpers.geometry.shapely_adapter import geometric_union


def phase_shifter_electrodes(
    cell,
    reference_port,
    crossing_goal_position,
    direction,
    electrode_wf_layer,
    param
):
    bounds = (np.inf, np.inf, -np.inf, -np.inf)
    wf_line_bounds = (np.inf, np.inf, -np.inf, -np.inf)

    if np.isclose(reference_port.angle, np.pi/2):
        crossing_angle = np.pi
        if direction == -1:
            reference_port.origin[1] = (reference_port.origin[1]
                                        + param['electrode_length'])
    elif np.isclose(reference_port.angle, -np.pi/2):
        crossing_angle = 0
        if direction == 1:
            reference_port.origin[1] = (reference_port.origin[1]
                                        - param['electrode_length'])

    left_inport = Port(
        (
            reference_port.origin[0]
            - (param['electrode_sep'] / 2.0 + param['electrode_width'] / 2),
            reference_port.origin[1]
        ),
        direction*np.pi/2.0, param['electrode_width']
    )

    signal_electrode_width = min(param['electrode_width'],
                                 param['wg_sep'] - param['electrode_sep'])
    signal_inport = Port((reference_port.origin[0] + param['wg_sep']/2.0,
                          reference_port.origin[1]),
                         direction*np.pi/2.0, signal_electrode_width)
    right_inport = Port(
        (
            reference_port.origin[0]
            + param['wg_sep'] + param['electrode_sep'] / 2.0
            + param['electrode_width'] / 2.0,
            reference_port.origin[1]
        ),
        direction*np.pi/2.0, param['electrode_width']
    )

    signal_wg = Waveguide.make_at_port(signal_inport)

    if crossing_angle == 0:
        long_wg = Waveguide.make_at_port(left_inport)
        short_wg = Waveguide.make_at_port(right_inport)
    else:
        long_wg = Waveguide.make_at_port(right_inport)
        short_wg = Waveguide.make_at_port(left_inport)

    # short waveguide
    # the shorter waveguide must be shorter
    # to accomodate for the other two waveguides
    short_wg.add_straight_segment(param['electrode_length']
                                  - 2 * param['electrode_sep_y']
                                  - 2 * param['electrode_width'])
    # 90 degree turn
    short_wg._current_port.angle = crossing_angle
    short_wg._current_port.origin[0] = (short_wg.current_port.origin[0]
                                        - np.cos(crossing_angle)
                                        * short_wg.current_port.width / 2.0)
    short_wg._current_port.width = param['crossing_width']
    short_wg._current_port.origin[1] = (short_wg.current_port.origin[1]
                                        - direction
                                        * short_wg.current_port.width/2.0)
    wf_line_bounds = update_bounds(wf_line_bounds,
                                   short_wg.get_shapely_outline().bounds)
    short_wg.add_straight_segment(abs(crossing_goal_position
                                      - short_wg.current_port.origin[0]))
    short_wg.add_straight_segment(param['electrode_taper_length'],
                                  param['electrode_width'])
    bounds = update_bounds(bounds, short_wg.get_shapely_outline().bounds)

    # signal waveguide (always in the middle)
    # the signal waveguide must also be slightly shorter
    # to accomodate for the longest waveguide
    signal_wg.add_straight_segment(param['electrode_length']
                                   - param['electrode_sep_y']
                                   - param['electrode_width'])
    # 90 degree turn
    signal_wg._current_port.angle = crossing_angle
    signal_wg._current_port.origin[0] = (signal_wg.current_port.origin[0]
                                         - np.cos(crossing_angle)
                                         * signal_wg.current_port.width/2.0)
    signal_wg._current_port.width = param['crossing_width']
    signal_wg._current_port.origin[1] = (signal_wg.current_port.origin[1]
                                         - direction
                                         * signal_wg.current_port.width/2.0)
    wf_line_bounds = update_bounds(wf_line_bounds,
                                   signal_wg.get_shapely_outline().bounds)
    # crossing segment + taper
    signal_wg.add_straight_segment(abs(crossing_goal_position
                                       - signal_wg.current_port.origin[0]))
    signal_wg.add_straight_segment(param['electrode_taper_length'],
                                   param['electrode_width'])
    bounds = update_bounds(bounds, signal_wg.get_shapely_outline().bounds)

    # long waveguide
    long_wg.add_straight_segment(param['electrode_length'])
    # 90 degree turn
    long_wg._current_port.angle = crossing_angle
    long_wg._current_port.origin[0] = (long_wg.current_port.origin[0]
                                       - np.cos(crossing_angle)
                                       * long_wg.current_port.width/2.0)
    long_wg._current_port.width = param['crossing_width']
    long_wg._current_port.origin[1] = (long_wg.current_port.origin[1]
                                       - direction
                                       * long_wg.current_port.width/2.0)
    wf_line_bounds = update_bounds(wf_line_bounds,
                                   long_wg.get_shapely_outline().bounds)
    # crossing segment + taper
    long_wg.add_straight_segment(abs(crossing_goal_position
                                     - long_wg.current_port.origin[0]))
    long_wg.add_straight_segment(param['electrode_taper_length'],
                                 param['electrode_width'])
    bounds = update_bounds(bounds, long_wg.get_shapely_outline().bounds)
    # the extra width ensures no tapers are clipped
    extra_width = (param['electrode_width']/2
                   - param['crossing_width']/2)
    wf_line_bounds = (wf_line_bounds[0],
                      wf_line_bounds[1] - extra_width,
                      wf_line_bounds[2],
                      wf_line_bounds[3] + extra_width)
    _, _ = wf_line_from_bounds(
        cell=cell,
        bounds=wf_line_bounds,
        wf_maxlength=param['wf_maxlength'],
        wf_layer=electrode_wf_layer
    )

    cell.add_to_layer(param['electrode_layer'], short_wg, signal_wg, long_wg)
    return (short_wg, signal_wg, long_wg, direction, wf_line_bounds)


def build_kinks(cell,
                left_wg,
                signal_wg,
                right_wg,
                kink,
                region_marker,
                electrode_wf_layer,
                previous_line_bounds,
                left_side,
                param):
    initial_angle = left_wg.current_port.angle
    # left_wg
    if left_side:
        kink_1_line_bounds = [np.inf,
                              np.inf,
                              previous_line_bounds[0],
                              -np.inf]
        wf_direction = -1
    else:
        kink_1_line_bounds = [previous_line_bounds[2],
                              np.inf,
                              -np.inf,
                              -np.inf]
        wf_direction = 1
    left_wg.add_straight_segment_until_x(kink[0]
                                         - param['electrode_pitch'])
    if left_side:
        kink_1_line_bounds[0] = (left_wg.current_port.origin[0]
                                 - param['wf_leeways'][0]/2)
        kink_1_line_bounds[1] = (left_wg.current_port.origin[1]
                                 - left_wg.current_port.width/2)
    else:
        kink_1_line_bounds[3] = (left_wg.current_port.origin[1]
                                 + left_wg.current_port.width/2)
    # 90 degree turn
    left_wg._current_port.origin[0] = (left_wg.current_port.origin[0]
                                       - np.cos(left_wg.current_port.angle)
                                       * left_wg.current_port.width/2.0)
    left_wg._current_port.angle = np.pi/2
    # straight segment to the next turn
    left_wg.add_straight_segment_until_y(kink[1]
                                         + np.cos(initial_angle)
                                         * param['electrode_pitch'])
    # 90 degree turn
    left_wg._current_port.angle = initial_angle
    left_wg._current_port.origin[0] = (left_wg._current_port.origin[0]
                                       - np.cos(initial_angle)
                                       * left_wg.current_port.width/2)

    # signal_wg
    signal_wg.add_straight_segment_until_x(kink[0])
    # 90 degree turn
    signal_wg._current_port.origin[0] = (signal_wg.current_port.origin[0]
                                         - np.cos(signal_wg.current_port.angle)
                                         * signal_wg.current_port.width/2.0)
    signal_wg._current_port.angle = np.pi/2
    # straight segment to the next turn
    signal_wg.add_straight_segment_until_y(kink[1])
    # 90 degree turn
    signal_wg._current_port.angle = initial_angle
    signal_wg._current_port.origin[0] = (signal_wg._current_port.origin[0]
                                         - np.cos(initial_angle)
                                         * signal_wg.current_port.width/2)

    # right_wg
    right_wg.add_straight_segment_until_x(kink[0]
                                          + param['electrode_pitch'])
    if not left_side:
        kink_1_line_bounds[2] = (right_wg.current_port.origin[0]
                                 + param['wf_leeways'][0]/2)
        kink_1_line_bounds[1] = (right_wg.current_port.origin[1]
                                 - right_wg.current_port.width/2)
    else:
        kink_1_line_bounds[3] = (right_wg.current_port.origin[1]
                                 + right_wg.current_port.width/2)
    # 90 degree turn
    right_wg._current_port.origin[0] = (right_wg.current_port.origin[0]
                                        - np.cos(right_wg.current_port.angle)
                                        * right_wg.current_port.width/2.0)
    right_wg._current_port.angle = np.pi/2
    # straight segment to the next turn
    right_wg.add_straight_segment_until_y(kink[1]
                                          - np.cos(initial_angle)
                                          * param['electrode_pitch'])
    # 90 degree turn
    right_wg._current_port.angle = initial_angle
    right_wg._current_port.origin[0] = (right_wg._current_port.origin[0]
                                        - np.cos(initial_angle)
                                        * right_wg.current_port.width/2)
    _, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=kink_1_line_bounds,
        region_marker=region_marker,
        wf_maxlength=param['wf_maxlength'],
        wf_layer=electrode_wf_layer,
        axis=0,
        direction=wf_direction
    )
    if left_side:
        kink_2_line_bounds = [
            kink_1_line_bounds[0],
            kink_1_line_bounds[3],
            (right_wg.current_port.origin[0]
             + param['wf_leeways'][0]/2),
            (left_wg.current_port.origin[1]
             - left_wg.current_port.width/2
             - param['wf_leeways'][1]/2)
        ]
        _, region_marker = wf_line_from_bounds(
            cell=cell,
            bounds=kink_2_line_bounds,
            region_marker=region_marker,
            wf_maxlength=param['wf_maxlength'],
            wf_layer=electrode_wf_layer
        )
        kink_2_line_bounds[0] = kink_2_line_bounds[2]
        # Tells where the next electrodes should branch in order
        # to run alongside this one
        next_kink = [kink[0] - 3 * param['electrode_pitch'],
                     left_wg.current_port.origin[1]
                     - 2 * param['electrode_pitch']]
    else:
        kink_2_line_bounds = [
            left_wg.current_port.origin[0] - param['wf_leeways'][0]/2,
            kink_1_line_bounds[3],
            (right_wg.current_port.origin[0]
             + right_wg.current_port.width
             + param['wf_leeways'][0]/2),
            (right_wg.current_port.origin[1]
             - right_wg.current_port.width/2
             - param['wf_leeways'][1]/2)
        ]
        _, region_marker = wf_line_from_bounds(
            cell=cell,
            bounds=kink_2_line_bounds,
            region_marker=region_marker,
            wf_maxlength=param['wf_maxlength'],
            wf_layer=electrode_wf_layer
        )
        kink_2_line_bounds[2] = kink_2_line_bounds[0]
        # Tells where the next electrodes should branch in order
        # to run alongside this one
        next_kink = [kink[0] + 3 * param['electrode_pitch'],
                     right_wg.current_port.origin[1]
                     - 2 * param['electrode_pitch']]
    return kink_2_line_bounds, next_kink, region_marker


def build_wrap_kinks(
    cell,
    left_wg,
    signal_wg,
    right_wg,
    kink,
    region_marker,
    wg_top,
    electrode_wf_layer,
    previous_line_bounds,
    left_side,
    param
):
    if kink == [-1, -1]:
        if left_side:
            kink = [(signal_wg.current_port.origin[0]
                     - 2 * param['electrode_pitch']),
                    (wg_top + 2*param['electrode_pitch'])]
        else:
            kink = [(signal_wg.current_port.origin[0]
                     + 2 * param['electrode_pitch']),
                    (wg_top + 2*param['electrode_pitch'])]
    initial_angle = left_wg.current_port.angle
    if left_side:
        kink_1_line_bounds = [np.inf,
                              np.inf,
                              previous_line_bounds[0],
                              -np.inf]
        wf_direction = -1
    else:
        kink_1_line_bounds = [previous_line_bounds[2],
                              np.inf,
                              -np.inf,
                              -np.inf]
        wf_direction = 1
    left_wg.add_straight_segment_until_x(kink[0]
                                         - param['electrode_pitch'])
    if left_side:
        kink_1_line_bounds[0] = (left_wg.current_port.origin[0]
                                 - param['wf_leeways'][0]/2)
        kink_1_line_bounds[1] = (left_wg.current_port.origin[1]
                                 - left_wg.current_port.width/2)
    else:
        kink_1_line_bounds[3] = (left_wg.current_port.origin[1]
                                 + left_wg.current_port.width/2)
    # 90 degree turn
    left_wg._current_port.origin[0] = (left_wg.current_port.origin[0]
                                       - np.cos(left_wg.current_port.angle)
                                       * left_wg.current_port.width/2.0)
    left_wg._current_port.angle = np.pi/2
    # straight segment to the next turn
    left_wg.add_straight_segment_until_y(kink[1]
                                         - np.cos(initial_angle)
                                         * param['electrode_pitch'])
    # 90 degree turn
    left_wg._current_port.angle = initial_angle + np.pi
    left_wg._current_port.origin[0] = (left_wg._current_port.origin[0]
                                       + np.cos(initial_angle)
                                       * left_wg.current_port.width/2)

    # signal_wg
    signal_wg.add_straight_segment_until_x(kink[0])
    # 90 degree turn
    signal_wg._current_port.origin[0] = (signal_wg.current_port.origin[0]
                                         - np.cos(signal_wg.current_port.angle)
                                         * signal_wg.current_port.width/2.0)
    signal_wg._current_port.angle = np.pi/2
    # straight segment to the next turn
    signal_wg.add_straight_segment_until_y(kink[1])
    # 90 degree turn
    signal_wg._current_port.angle = initial_angle + np.pi
    signal_wg._current_port.origin[0] = (signal_wg._current_port.origin[0]
                                         + np.cos(initial_angle)
                                         * signal_wg.current_port.width/2)

    # right_wg
    right_wg.add_straight_segment_until_x(kink[0]
                                          + param['electrode_pitch'])
    if not left_side:
        kink_1_line_bounds[2] = (right_wg.current_port.origin[0]
                                 + param['wf_leeways'][0]/2)
        kink_1_line_bounds[1] = (right_wg.current_port.origin[1]
                                 - right_wg.current_port.width/2)
    else:
        kink_1_line_bounds[3] = (right_wg.current_port.origin[1]
                                 + right_wg.current_port.width/2)
    # 90 degree turn
    right_wg._current_port.origin[0] = (right_wg.current_port.origin[0]
                                        - np.cos(right_wg.current_port.angle)
                                        * right_wg.current_port.width/2.0)
    right_wg._current_port.angle = np.pi/2
    # straight segment to the next turn
    right_wg.add_straight_segment_until_y(kink[1]
                                          + np.cos(initial_angle)
                                          * param['electrode_pitch'])
    # 90 degree turn
    right_wg._current_port.angle = initial_angle + np.pi
    right_wg._current_port.origin[0] = (right_wg._current_port.origin[0]
                                        + np.cos(initial_angle)
                                        * right_wg.current_port.width/2)
    _, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=kink_1_line_bounds,
        region_marker=region_marker,
        wf_maxlength=param['wf_maxlength'],
        wf_layer=electrode_wf_layer,
        axis=0,
        direction=wf_direction
    )
    if not left_side:
        kink_2_line_bounds = [
            (left_wg.current_port.origin[0]
             - left_wg.current_port.width
             - param['wf_leeways'][0]/2),
            kink_1_line_bounds[3],
            (right_wg.current_port.origin[0]
             + param['wf_leeways'][0]/2),
            (left_wg.current_port.origin[1]
             - left_wg.current_port.width/2
             - param['wf_leeways'][1]/2)
        ]
        _, region_marker = wf_line_from_bounds(
            cell=cell,
            bounds=kink_2_line_bounds,
            region_marker=region_marker,
            wf_maxlength=param['wf_maxlength'],
            wf_layer=electrode_wf_layer
        )
        kink_2_line_bounds[0] = kink_2_line_bounds[2]
        # Tells where the next electrodes should branch in order
        # to run alongside this one
        next_kink = [kink[0] + 3 * param['electrode_pitch'],
                     right_wg.current_port.origin[1]
                     + 2 * param['electrode_pitch']]
    else:
        kink_2_line_bounds = [
            left_wg.current_port.origin[0] - param['wf_leeways'][0]/2,
            kink_1_line_bounds[3],
            (right_wg.current_port.origin[0]
             + right_wg.current_port.width
             + param['wf_leeways'][0]/2),
            (right_wg.current_port.origin[1]
             - right_wg.current_port.width/2
             - param['wf_leeways'][1]/2)
        ]
        _, region_marker = wf_line_from_bounds(
            cell=cell,
            bounds=kink_2_line_bounds,
            region_marker=region_marker,
            wf_maxlength=param['wf_maxlength'],
            wf_layer=electrode_wf_layer
        )
        kink_2_line_bounds[2] = kink_2_line_bounds[0]
        # Tells where the next electrodes should branch in order
        # to run alongside this one
        next_kink = [kink[0] - 3 * param['electrode_pitch'],
                     left_wg.current_port.origin[1]
                     + 2 * param['electrode_pitch']]
    return kink_2_line_bounds, next_kink, region_marker


def electrode_probe_and_pad(cell,
                            left_wg,
                            signal_wg,
                            right_wg,
                            region_marker,
                            electrode_wf_layer,
                            connector_coordinates,
                            ground_signal_pad_sep,
                            first_line_bounds,
                            left_side,
                            param):
    left_wg.add_straight_segment_until_y(
        connector_coordinates[1],
        final_width=param['connector_probe_dims'][0]
    )

    signal_wg.add_straight_segment_until_y(
        connector_coordinates[1],
        final_width=param['connector_probe_dims'][0]
    )
    right_wg.add_straight_segment_until_y(
        connector_coordinates[1],
        final_width=param['crossing_width']
    )
    # probe connector
    left_wg.add_straight_segment(
        param['connector_probe_dims'][1]
    )
    signal_wg.add_straight_segment(
        param['connector_probe_dims'][1]
    )
    right_wg.add_straight_segment(
        param['connector_probe_dims'][1]
    )
    if left_side:
        second_line_bounds = [(left_wg.current_port.origin[0]
                               - left_wg.current_port.width/2
                               - param['wf_leeways'][0]/2),
                              first_line_bounds[1],
                              first_line_bounds[0],
                              (left_wg.current_port.origin[1]
                               - param['electrode_width']
                               - param['wf_leeways'][1]/2)]
        if second_line_bounds[2] < (right_wg.current_port.origin[0]
                                    + right_wg.current_port.width/2):
            second_line_bounds[2] = (right_wg.current_port.origin[0]
                                     + right_wg.current_port.width/2
                                     + param['wf_leeways'][0]/2)
    else:
        second_line_bounds = [first_line_bounds[2],
                              first_line_bounds[1],
                              (right_wg.current_port.origin[0]
                               + param['electrode_width']/2
                               + param['wf_leeways'][0]/2),
                              (left_wg.current_port.origin[1]
                               - param['electrode_width']
                               - param['wf_leeways'][1]/2)]
        if second_line_bounds[0] > (left_wg.current_port.origin[0]
                                    - left_wg.current_port.width/2):
            second_line_bounds[0] = (left_wg.current_port.origin[0]
                                     - left_wg.current_port.width/2
                                     - param['wf_leeways'][0]/2)
    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=second_line_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
    )
    first_pad_wf_bounds = [-1,
                           second_line_bounds[3],
                           -1,
                           -1]
    second_pad_wf_bounds = [-1, -1, -1, -1]
    intermediate_pad_wf_bounds = [-1, -1, -1, -1]

    # electrode_connector
    # make sure the contact pad does not touch the ground electrodes
    signal_wg.add_straight_segment(param['electrode_pad_seps'][1])
    if (
        signal_wg.current_port.width
        > (param['contact_pad_dims'][0]/2
           - param['connector_probe_pitch'])
    ):
        pad_position = (
            signal_wg.current_port.origin[0]
            - param['connector_probe_pitch']/2,
            signal_wg.current_port.origin[1]
        )
        signal_wg.add_straight_segment(param['contact_pad_dims'][1])
        signal_wg._current_port.origin = pad_position
    else:
        signal_wg._current_port.origin[0] = (
            signal_wg.current_port.origin[0]
            - param['connector_probe_pitch']/2
        )
    signal_wg._current_port.width = param['contact_pad_dims'][0]
    signal_wg.add_straight_segment(param['contact_pad_dims'][1])

    # Route the left wg to the ground pad
    left_wg._current_port.width = param['crossing_width']
    left_probe_x = left_wg.current_port.origin[0]

    x_goal_length = abs((signal_wg.current_port.origin[0]
                         - param['contact_pad_dims'][0]/2
                         - ground_signal_pad_sep)
                        - left_wg.current_port.origin[0])
    if x_goal_length > param['connector_probe_dims'][0]/2:
        left_wg._current_port.angle = np.pi
        left_wg._current_port.origin[1] = (left_wg.current_port.origin[1]
                                           - left_wg.current_port.width/2)
        left_wg.add_straight_segment(x_goal_length)
        left_wg._current_port.angle = np.pi/2
        left_wg._current_port.origin[1] = (left_wg.current_port.origin[1]
                                           - left_wg.current_port.width/2)
    else:
        left_wg._current_port.origin[0] = (
            left_wg.current_port.origin[0]
            - param['connector_probe_dims'][0]/2
            + left_wg.current_port.width/2
        )

    left_wg.add_straight_segment_until_y(signal_wg.current_port.origin[1]
                                         + param['contact_pad_sep']
                                         + left_wg.current_port.width/2)
    first_pad_wf_bounds[0] = (left_wg.current_port.origin[0]
                              - left_wg.current_port.width/2
                              - param['wf_leeways'][0]/8)
    left_wg._current_port.origin[0] = (
        left_probe_x
        - param['connector_probe_pitch']/2
    )
    left_wg._current_port.width = param['contact_pad_dims'][0]
    left_wg.add_straight_segment(param['contact_pad_dims'][1])
    second_pad_wf_bounds[0] = (left_wg.current_port.origin[0]
                               - left_wg.current_port.width/2
                               - param['wf_leeways'][0])
    second_pad_wf_bounds[2] = (left_wg.current_port.origin[0]
                               + left_wg.current_port.width/2
                               + 38)
    intermediate_pad_wf_bounds[0] = first_pad_wf_bounds[0]
    intermediate_pad_wf_bounds[2] = second_pad_wf_bounds[2]

    right_wg.add_straight_segment_until_y(
        signal_wg.current_port.origin[1]
        + param['electrode_pad_seps'][0]
        + right_wg.current_port.width/2
    )
    first_pad_wf_bounds[2] = (right_wg.current_port.origin[0]
                              + right_wg.current_port.width/2
                              + param['wf_leeways'][0]/8)

    right_wg._current_port.angle = np.pi
    right_wg._current_port.origin[0] = (right_wg.current_port.origin[0]
                                        + right_wg.current_port.width/2)
    right_wg.add_straight_segment_until_x(left_wg.current_port.origin[0]
                                          + param['contact_pad_dims'][0]/2
                                          - right_wg.current_port.width/2)
    first_pad_wf_bounds[3] = (
        right_wg.current_port.origin[1]
        + right_wg.current_port.width
    )
    right_wg._current_port.angle = np.pi/2
    right_wg._current_port.origin[1] = (
        right_wg.current_port.origin[1]
        - right_wg.current_port.width/2
    )
    right_wg.add_straight_segment_until_y(signal_wg.current_port.origin[1]
                                          + param['contact_pad_sep']
                                          + right_wg.current_port.width/2)
    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=first_pad_wf_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
            )

    intermediate_pad_wf_bounds[1] = first_pad_wf_bounds[3]
    intermediate_pad_wf_bounds[3] = first_pad_wf_bounds[3] + 10
    second_pad_wf_bounds[1] = intermediate_pad_wf_bounds[3]
    second_pad_wf_bounds[3] = (first_pad_wf_bounds[3]
                               + 400)
    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=intermediate_pad_wf_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
            )
    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=second_pad_wf_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
            )
    cell.add_to_layer(electrode_wf_layer+1,
                      region_marker)


def electrode_pad(cell,
                  left_wg,
                  signal_wg,
                  right_wg,
                  region_marker,
                  electrode_wf_layer,
                  connector_coordinates,
                  first_line_bounds,
                  left_side,
                  param):
    left_wg.add_straight_segment_until_y(
        connector_coordinates[1]
    )
    signal_wg.add_straight_segment_until_y(
        connector_coordinates[1]
    )
    right_wg.add_straight_segment_until_y(
        connector_coordinates[1]
    )

    if left_side:
        second_line_bounds = [(left_wg.current_port.origin[0]
                               - left_wg.current_port.width/2
                               - param['wf_leeways'][0]/2),
                              first_line_bounds[1],
                              first_line_bounds[0],
                              (left_wg.current_port.origin[1]
                               - param['electrode_width']
                               - param['wf_leeways'][1]/2)]
        if second_line_bounds[2] < (right_wg.current_port.origin[0]
                                    + right_wg.current_port.width/2):
            second_line_bounds[2] = (right_wg.current_port.origin[0]
                                     + right_wg.current_port.width/2
                                     + param['wf_leeways'][0]/2)
    else:
        second_line_bounds = [first_line_bounds[2],
                              first_line_bounds[1],
                              (right_wg.current_port.origin[0]
                               + right_wg.current_port.width/2
                               + param['wf_leeways'][0]/2),
                              (left_wg.current_port.origin[1]
                               - param['electrode_width']
                               - param['wf_leeways'][1]/2)]
        if second_line_bounds[0] > (left_wg.current_port.origin[0]
                                    - left_wg.current_port.width/2):
            second_line_bounds[0] = (left_wg.current_port.origin[0]
                                     - left_wg.current_port.width/2
                                     - param['wf_leeways'][0]/2)

    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=second_line_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
            )
    first_pad_wf_bounds = [-1,
                           second_line_bounds[3],
                           -1,
                           -1]
    second_pad_wf_bounds = [-1, -1, -1, -1]

    # electrode_connector
    # make sure the contact pad does not touch the ground electrodes
    signal_wg.add_straight_segment(param['electrode_pad_seps'][1])
    signal_wg._current_port.width = param['contact_pad_dims'][0]
    signal_wg.add_straight_segment(param['contact_pad_dims'][1])

    left_wg.add_straight_segment_until_y(signal_wg.current_port.origin[1]
                                         + param['contact_pad_sep']
                                         + left_wg.current_port.width/2)
    first_pad_wf_bounds[0] = (left_wg.current_port.origin[0]
                              - left_wg.current_port.width/2
                              - param['wf_leeways'][0])
    left_wg._current_port.width = param['contact_pad_dims'][0]
    left_wg.add_straight_segment(param['contact_pad_dims'][1])
    second_pad_wf_bounds[0] = (left_wg.current_port.origin[0]
                               - left_wg.current_port.width/2
                               - param['wf_leeways'][0])
    second_pad_wf_bounds[2] = (left_wg.current_port.origin[0]
                               + left_wg.current_port.width/2
                               + 60)
    right_wg._current_port.width = param['electrode_width']

    right_wg.add_straight_segment_until_y(
        signal_wg.current_port.origin[1]
        + param['electrode_pad_seps'][0]
        + right_wg.current_port.width/2
    )
    first_pad_wf_bounds[2] = (right_wg.current_port.origin[0]
                              + right_wg.current_port.width/2
                              + param['wf_leeways'][0])

    right_wg._current_port.angle = np.pi
    right_wg._current_port.origin[0] = (right_wg.current_port.origin[0]
                                        + right_wg.current_port.width/2)
    right_wg.add_straight_segment_until_x(left_wg.current_port.origin[0]
                                          + param['contact_pad_dims'][0]/2
                                          - right_wg.current_port.width/2)
    first_pad_wf_bounds[3] = (
        right_wg.current_port.origin[1]
        + right_wg.current_port.width
    )
    right_wg._current_port.angle = np.pi/2
    right_wg._current_port.origin[1] = (
        right_wg.current_port.origin[1]
        - right_wg.current_port.width/2
    )
    right_wg.add_straight_segment_until_y(signal_wg.current_port.origin[1]
                                          + param['contact_pad_sep']
                                          + right_wg.current_port.width/2)
    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=first_pad_wf_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
            )
    second_pad_wf_bounds[1] = first_pad_wf_bounds[3]
    second_pad_wf_bounds[3] = (first_pad_wf_bounds[3]
                               + 650)
    _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=second_pad_wf_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer
            )
    # add markers
    cell.add_to_layer(electrode_wf_layer+1,
                      region_marker)


def electrode_connector(cell,
                        left_wg,
                        signal_wg,
                        right_wg,
                        kink,
                        wg_top,
                        electrode_wf_layer,
                        region_marker,
                        connector_coordinates,
                        previous_line_bounds,
                        left_side,
                        initial_point,
                        param,
                        probe_connector=True):
    # Collect electrodes to save space
    next_kink = [-1, -1]
    new_region_marker = bounds_to_polygon(previous_line_bounds)
    if left_side:
        end_point = (previous_line_bounds[2],
                     (previous_line_bounds[1] + previous_line_bounds[3]) / 2)
    else:
        end_point = (previous_line_bounds[0],
                     (previous_line_bounds[1] + previous_line_bounds[3]) / 2)
    region_marker = connect_writefields(cell,
                                        initial_point,
                                        end_point,
                                        electrode_wf_layer,
                                        region_marker,
                                        param)
    region_marker = geometric_union([
            region_marker,
            new_region_marker
    ])

    if kink != [-1, -1]:
        if (left_side and not (
                kink[0] + param['electrode_pitch']
                > (connector_coordinates[0] + param['connector_probe_pitch'])
        )):
            previous_line_bounds, next_kink, region_marker = build_wrap_kinks(
                cell,
                left_wg,
                signal_wg,
                right_wg,
                kink,
                region_marker,
                wg_top,
                electrode_wf_layer,
                previous_line_bounds,
                left_side,
                param
            )
            left_side = False
        elif (not left_side
              and not (
                  kink[0] - param['electrode_pitch'] <
                  connector_coordinates[0] - param['connector_probe_pitch']
                  )
              ):
            previous_line_bounds, next_kink, region_marker = build_wrap_kinks(
                cell,
                left_wg,
                signal_wg,
                right_wg,
                kink,
                region_marker,
                wg_top,
                electrode_wf_layer,
                previous_line_bounds,
                left_side,
                param
            )
            left_side = True
        else:
            previous_line_bounds, next_kink, region_marker = build_kinks(
                cell,
                left_wg,
                signal_wg,
                right_wg,
                kink,
                region_marker,
                electrode_wf_layer,
                previous_line_bounds,
                left_side,
                param
            )
    else:
        if (
            left_side
            and connector_coordinates[0] > signal_wg.current_port.origin[0]
        ):
            previous_line_bounds, next_kink, region_marker = build_wrap_kinks(
                cell,
                left_wg,
                signal_wg,
                right_wg,
                kink,
                region_marker,
                wg_top,
                electrode_wf_layer,
                previous_line_bounds,
                left_side,
                param
            )
            left_side = False
        elif (
            not left_side
            and connector_coordinates[0] < signal_wg.current_port.origin[0]
        ):
            previous_line_bounds, next_kink, region_marker = build_wrap_kinks(
                cell,
                left_wg,
                signal_wg,
                right_wg,
                kink,
                region_marker,
                wg_top,
                electrode_wf_layer,
                previous_line_bounds,
                left_side,
                param
            )
            left_side = True

    if left_side:
        first_line_bounds = [-1,
                             (left_wg.current_port.origin[1]
                              - left_wg.current_port.width/2
                              - param['wf_leeways'][1]/2),
                             previous_line_bounds[0],
                             (right_wg.current_port.origin[1]
                              + right_wg.current_port.width/2
                              + param['wf_leeways'][1]/2)]
        wf_direction = -1
    else:
        first_line_bounds = [previous_line_bounds[2],
                             (right_wg.current_port.origin[1]
                              - right_wg.current_port.width/2
                              - param['wf_leeways'][1]/2),
                             -1,
                             (left_wg.current_port.origin[1]
                              + left_wg.current_port.width/2
                              + param['wf_leeways'][1]/2)]
        wf_direction = 1

    # left_wg
    if probe_connector:
        ground_signal_pad_sep = 0.5 * 3/5 * (
            2 * param['connector_probe_pitch']
            - param['contact_pad_dims'][0]
        )
        left_wg_goal_x = (connector_coordinates[0]
                          - param['connector_probe_pitch'])
    else:
        left_wg_goal_x = (connector_coordinates[0]
                          - param['contact_pad_dims'][0]/2
                          - param['electrode_pad_seps'][0]
                          - left_wg.current_port.width/2)
    left_wg.add_straight_segment_until_x(left_wg_goal_x)
    if not left_side:
        first_line_bounds[2] = (left_wg.current_port.origin[0]
                                - left_wg.current_port.width/2
                                - param['connector_probe_dims'][0]/2
                                - param['wf_leeways'][0]/8)
        if (
            first_line_bounds[2] - first_line_bounds[0]
            < param['electrode_pitch'] + param['electrode_width']
           ):
            first_line_bounds[2] = first_line_bounds[0]
        else:
            _, region_marker = wf_line_from_bounds(
                    cell=cell,
                    bounds=first_line_bounds,
                    region_marker=region_marker,
                    wf_maxlength=param['wf_maxlength'],
                    wf_layer=electrode_wf_layer,
                    axis=0,
                    direction=wf_direction
            )
    # 90 degree turn
    left_wg._current_port.angle = np.pi/2
    left_wg._current_port.origin[1] = (left_wg._current_port.origin[1]
                                       - left_wg.current_port.width/2.0)

    # signal_wg (which is also the center wg)
    signal_wg_goal_x = connector_coordinates[0]
    signal_wg.add_straight_segment_until_x(signal_wg_goal_x)
    # 90 degree turn
    signal_wg._current_port.angle = np.pi/2
    signal_wg._current_port.origin[1] = (signal_wg._current_port.origin[1]
                                         - signal_wg.current_port.width/2.0)
    # right_wg
    if probe_connector:
        right_wg_goal_x = (connector_coordinates[0]
                           - param['connector_probe_pitch']/2
                           + param['contact_pad_dims'][0]/2
                           + ground_signal_pad_sep)
    else:
        right_wg_goal_x = (connector_coordinates[0]
                           + param['contact_pad_dims'][0]/2
                           + param['electrode_pad_seps'][0]
                           + right_wg.current_port.width/2)
    right_wg.add_straight_segment_until_x(right_wg_goal_x)
    if left_side:
        first_line_bounds[0] = (right_wg.current_port.origin[0]
                                + right_wg.current_port.width/2
                                + param['wf_leeways'][0]/8)
        if (
            first_line_bounds[2] - first_line_bounds[0]
            < param['electrode_pitch'] + param['electrode_width']
           ):
            first_line_bounds[0] = first_line_bounds[2]
            first_line_bounds[0] += param['electrode_width']/2
        else:
            _, region_marker = wf_line_from_bounds(
                cell=cell,
                bounds=first_line_bounds,
                region_marker=region_marker,
                wf_maxlength=param['wf_maxlength'],
                wf_layer=electrode_wf_layer,
                axis=0,
                direction=wf_direction
            )
    # 90 degree turn
    right_wg._current_port.origin[1] = (right_wg._current_port.origin[1]
                                        - right_wg.current_port.width/2.0)
    right_wg._current_port.angle = np.pi/2
    # Tells where the next electrodes should branch in order
    # to run alongside this one
    if next_kink == [-1, -1]:
        if left_side:
            next_kink = [right_wg_goal_x - param['electrode_pitch'],
                         left_wg.current_port.origin[1]
                         - 2 * param['electrode_pitch']
                         - left_wg.current_port.width/2.0]
        else:
            next_kink = [left_wg_goal_x + param['electrode_pitch'],
                         right_wg.current_port.origin[1]
                         - 2 * param['electrode_pitch']]
    if probe_connector:
        electrode_probe_and_pad(
            cell,
            left_wg,
            signal_wg,
            right_wg,
            region_marker,
            electrode_wf_layer,
            connector_coordinates,
            ground_signal_pad_sep,
            first_line_bounds,
            left_side,
            param
        )
    else:
        electrode_pad(
            cell,
            left_wg,
            signal_wg,
            right_wg,
            region_marker,
            electrode_wf_layer,
            connector_coordinates,
            first_line_bounds,
            left_side,
            param
        )
    return next_kink
