import numpy as np
from utils import *
from gdshelpers.parts.waveguide import Waveguide
from gdshelpers.parts.port import Port
from gdshelpers.parts.coupler import GratingCoupler
from euler_curves import wgAdd_EulerWiggle
from electrode_parts import electrode_connector, phase_shifter_electrodes
from wg_parts import wgs_to_fiber_array, expand_wgs_section_2
from wg_parts import phase_shifter_and_dc_wg

# Version 2.1

def marker_L_positions(center, dir_1, dir_2, separation=200):
    positions = [(center[0] + idx * dir_1 * separation, center[1])
                 for idx in range(3)]
    positions_ax_2 = [(center[0], center[1] + idx*dir_2*separation)
                      for idx in range(3)]
    positions[len(positions):] = positions_ax_2
    return positions


def put_markers_at_positions(cell, positions, marker_diameter, layer):
    markers = [SquareMarker.make_marker(position, marker_diameter)
               for position in positions]
    for marker in markers:
        cell.add_to_layer(layer, marker)


def build_global_markers(cell,
                         bounds,
                         marker_diameter,
                         marker_layer,
                         wg_layer,
                         marker_protection_layer):
    top_left = marker_L_positions(
        center=(bounds[0], bounds[3]),
        dir_1=1,
        dir_2=-1
    )
    put_markers_at_positions(cell, top_left, marker_diameter, marker_layer)
    put_markers_at_positions(cell, top_left, marker_diameter*1.5, wg_layer)
    put_markers_at_positions(cell, top_left,
                             marker_diameter*2, marker_protection_layer)
    top_right = marker_L_positions(
        center=(bounds[2], bounds[3]),
        dir_1=-1,
        dir_2=-1
    )
    put_markers_at_positions(cell, top_right, marker_diameter, marker_layer)
    put_markers_at_positions(cell, top_right, marker_diameter*1.5, wg_layer)
    put_markers_at_positions(cell, top_right,
                             marker_diameter*2, marker_protection_layer)
    bottom_right = marker_L_positions(
        center=(bounds[2], bounds[1]),
        dir_1=-1,
        dir_2=1
    )
    put_markers_at_positions(cell, bottom_right, marker_diameter, marker_layer)
    put_markers_at_positions(cell, bottom_right, marker_diameter*1.5, wg_layer)
    put_markers_at_positions(cell, bottom_right,
                             marker_diameter*2, marker_protection_layer)
    bottom_left = marker_L_positions(
        center=(bounds[0], bounds[1]),
        dir_1=1,
        dir_2=1
    )
    put_markers_at_positions(cell, bottom_left, marker_diameter, marker_layer)
    put_markers_at_positions(cell, bottom_left, marker_diameter*1.5, wg_layer)
    put_markers_at_positions(cell, bottom_left,
                             marker_diameter*2, marker_protection_layer)


def section_1_dcs_electrodes(
    cell,
    wgs,
    electrode_wf_layers,
    ending_taper,
    param
):
    # First two waveguides
    if wgs[0].current_port.angle == np.pi/2:
        direction_1 = -1
        direction_2 = 1
        crossing_goal_position = (wgs[0].current_port.origin[0]
                                  - wgs[0].current_port.width/2)
    else:
        direction_1 = 1
        direction_2 = -1
        crossing_goal_position = (wgs[-1].current_port.origin[0]
                                  + wgs[-1].current_port.width/2)

    stagger_separation = (2*(param['sine_s_x'] + param['coupler_length'])
                          + 3*param['electrode_width']
                          + 2*param['electrode_sep_y'])/2

    electrodes_1 = phase_shifter_electrodes(
        cell=cell,
        reference_port=wgs[0].current_port,
        crossing_goal_position=crossing_goal_position,
        direction=direction_1,
        electrode_wf_layer=electrode_wf_layers[0],
        param=param
    )
    phase_shifter_and_dc_wg(wgs[:2],
                            stagger_separation=stagger_separation,
                            ending_taper=ending_taper,
                            param=param)

    # Second two waveguides
    for wg in wgs[2:]:
        wg.add_straight_segment(stagger_separation)

    electrodes_2 = phase_shifter_electrodes(
        cell=cell,
        reference_port=wgs[2].current_port,
        crossing_goal_position=crossing_goal_position,
        direction=direction_2,
        electrode_wf_layer=electrode_wf_layers[1],
        param=param)
    phase_shifter_and_dc_wg(wgs[2:],
                            stagger_separation=0,
                            ending_taper=ending_taper,
                            param=param)
    return [electrodes_1, electrodes_2]


def section_2_dcs_electrodes(
    cell,
    wgs,
    electrode_wf_layer,
    ending_taper,
    param,
    delay=True,
    extra_wiggle_lengths=[-1, -1, -1, -1]
):
    # Interfering waveguides
    if wgs[1].current_port.angle == np.pi/2:
        direction_1 = 1
        crossing_goal_position = (wgs[0].current_port.origin[0]
                                  - wgs[0].current_port.width/2)
    else:
        direction_1 = -1
        crossing_goal_position = (wgs[-1].current_port.origin[0]
                                  + wgs[-1].current_port.width/2)
    electrodes = phase_shifter_electrodes(
        cell=cell,
        reference_port=wgs[1].current_port,
        crossing_goal_position=crossing_goal_position,
        direction=direction_1,
        electrode_wf_layer=electrode_wf_layer,
        param=param)
    phase_shifter_and_dc_wg(wgs[1:3],
                            ending_taper=ending_taper,
                            param=param)

    # Delayed waveguides
    if not delay:
        return electrodes
    for idx in (0, 3):
        current_y_position = wgs[idx].current_port.origin[1]
        goal_y_position = wgs[1].current_port.origin[1]

        if extra_wiggle_lengths == [-1, -1, -1, -1]:
            path_length_difference = wgs[1].length - wgs[idx].length
        else:
            path_length_difference = (
                extra_wiggle_lengths[idx]
                + abs(goal_y_position-current_y_position)
            )

        radius = param['min_radius']
        target_path_length = path_length_difference
        target_crow_length = abs(goal_y_position-current_y_position)
        if target_path_length > target_crow_length:
            wgAdd_EulerWiggle(wgs[idx],
                              radius=radius,
                              target_path_length=target_path_length,
                              target_crow_length=target_crow_length,
                              internal_angle_mod=0.0,
                              N_turns=3,
                              mirrored=False,
                              resolution=200)
        else:
            wgs[idx].add_straight_segment(
                abs(goal_y_position - current_y_position)
            )

    return electrodes


def build_interferometer(
    cell,
    wgs,
    right_electrode_region_markers,
    left_electrode_region_markers,
    region_marker,
    second_x_middle,
    initial_wf_point,
    delay,
    param,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf)
):
    initial_position_y = wgs[0].current_port.origin[1]
    x_middle = (wgs[0].current_port.origin[0]
                + wgs[-1].current_port.origin[0])/2
    wg_wf_bounds = list(total_bounds).copy()

    # connector probe centers start coordinates
    global_x_middle = (x_middle + second_x_middle)/2
    connector_center_separation = 2 * param['connector_probe_pitch']
    # Center slightly to the left as there are more electrodes to the right
    left_connector_x_middle = (global_x_middle
                               - connector_center_separation)
    right_connector_x_middle = (global_x_middle)
    left_connector_coordinates = [(left_connector_x_middle
                                   - (3-idx)*connector_center_separation,
                                   param['connector_start_y'])
                                  for idx in range(4)]
    right_connector_coordinates = [(right_connector_x_middle
                                    + idx*connector_center_separation,
                                    param['connector_start_y'])
                                   for idx in range(6)]

    # assume the wg separation is already appropriate for section 1
    # First section 1
    left_side_electrodes = [None for _ in range(5)]
    right_side_electrodes = [None for _ in range(6)]
    # straight segment for the electrode with taper
    stagger_separation = (2*(param['sine_s_x'] + param['coupler_length'])
                          + 3*param['electrode_width']
                          + 2*param['electrode_sep_y'])/2
    for wg in wgs:
        wg.add_straight_segment(param['mm_taper_length'],
                                param['mm_wg_width'])
        wg.add_straight_segment(param['electrode_length']
                                + stagger_separation
                                + param['coupler_length']
                                + 2*param['sine_s_x']
                                + param['mm_taper_length'])

    for electrode_idx in range(1):
        sec1_electrodes = section_1_dcs_electrodes(
            cell=cell,
            wgs=wgs,
            electrode_wf_layers=param['left_electrode_wf_layers'][
                2*electrode_idx:2*(electrode_idx+2)
            ],
            ending_taper=(electrode_idx) % 2,
            param=param
        )
        left_side_electrodes[2*electrode_idx] = sec1_electrodes[0]
        left_side_electrodes[2*electrode_idx+1] = sec1_electrodes[1]
    # First section 2
    # Expand the wgs first
    _, _ = expand_wgs_section_2(wgs, param)
    for electrode_idx in range(2, 4):
        sec2_electrodes = section_2_dcs_electrodes(
            cell=cell,
            wgs=wgs,
            electrode_wf_layer=(
                param['left_electrode_wf_layers'][electrode_idx]
            ),
            ending_taper=(electrode_idx+1) % 2,
            param=param,
            delay=False
        )
        left_side_electrodes[electrode_idx] = sec2_electrodes

    # first and last mode delay lines
    curve_goal_x_positions = [
        # leftmost to rightmost
        (second_x_middle
         + 0.5*param['electrode_wg_sep']
         + param['wg_sep']),
        # #2 from left to #2 from right
        (second_x_middle
         + 0.5*param['electrode_wg_sep']),
        # #2 from right to #2 from left
        (second_x_middle
         - 0.5*param['electrode_wg_sep']),
        # rightmost to leftmost
        (second_x_middle
         - 0.5*param['electrode_wg_sep']
         - param['wg_sep'])
    ]
    wg_sep = 90
    left_to_right = wgs[0]
    if param['mm_taper_length'] > 0:
        left_to_right.add_straight_segment(param['mm_taper_length'],
                                           param['mm_wg_width'])
    left_to_right.add_straight_segment_until_y(
        wgs[1].current_port.origin[1]
        + 3 * wg_sep
        - 50
        - param['mm_taper_length']
    )
    if param['mm_taper_length'] > 0:
        left_to_right.add_straight_segment(param['mm_taper_length'],
                                           param['sm_wg_width'])
    left_to_right.add_bend(-np.pi/2, param['min_radius'])
    if param['mm_taper_length'] > 0:
        left_to_right.add_straight_segment(param['mm_taper_length'],
                                           param['mm_wg_width'])
    left_to_right.add_straight_segment_until_x(
        curve_goal_x_positions[0]
        - param['min_radius']
        - param['mm_taper_length']
    )
    if param['mm_taper_length'] > 0:
        left_to_right.add_straight_segment(param['mm_taper_length'],
                                           param['sm_wg_width'])
    left_to_right.add_bend(-np.pi/2, param['min_radius'])
    left_to_right.add_straight_segment(3*wg_sep - 50)

    longest_path_length = left_to_right.length
    path_length = [0, 0, 0]
    initial_crow_length = (
        wgs[1].current_port.origin[1]
        - wgs[3].current_port.origin[1]
    )

    for idx in range(1, 4):
        curve_length = abs(
            curve_goal_x_positions[idx]
            - wgs[idx].current_port.origin[0]
        )
        path_length[idx-1] = (
            wgs[idx].length
            + np.pi * curve_length
        )
        if idx == 3:
            if delay:
                radius = param['min_radius']
                # initial_length = wgs[3].length
                target_path_length = (
                    longest_path_length
                    - wgs[3].length
                    - np.pi * param['min_radius']
                    - curve_goal_x_positions[3]
                    + wgs[3].current_port.origin[0]
                    + 2*param['min_radius']
                )
                target_crow_length = initial_crow_length
                wgAdd_EulerWiggle(
                    wgs[3],
                    radius=radius,
                    target_path_length=target_path_length,
                    target_crow_length=target_crow_length,
                    internal_angle_mod=0.0,
                    N_turns=49,
                    resolution=200,
                    mirrored=True
                )
                first_line_right_bound = (
                    wgs[3].get_shapely_outline().bounds[2]
                )
            else:
                first_line_right_bound = (wgs[3].current_port.origin[0]
                                          + wgs[3].current_port.width)
                wgs[3].add_straight_segment(initial_crow_length)
            # final_length = wgs[3].length
            wgs[3].add_bend(-np.pi/2, param['min_radius'])
            if param['mm_taper_length'] > 0:
                wgs[3].add_straight_segment(param['mm_taper_length'],
                                            param['mm_wg_width'])
            wgs[3].add_straight_segment_until_x(
                curve_goal_x_positions[3]
                - param['min_radius']
                - param['mm_taper_length']
            )
            if param['mm_taper_length'] > 0:
                wgs[3].add_straight_segment(param['mm_taper_length'],
                                            param['sm_wg_width'])
            wgs[3].add_bend(-np.pi/2, param['min_radius'])

        else:
            wgs[idx].add_straight_segment((3-idx) * wg_sep - 50)
            wgs[idx].add_bend(-np.pi/2, param['min_radius'])
            target_path_length = (longest_path_length
                                  - wgs[idx].length
                                  - np.pi/2 * param['min_radius']
                                  - (3-idx) * wg_sep
                                  + 50)
            target_crow_length = (curve_goal_x_positions[idx]
                                  - wgs[idx].current_port.origin[0]
                                  - param['min_radius'])
            if delay:
                radius = param['min_radius']
                wgAdd_EulerWiggle(
                    wgs[idx],
                    radius=radius,
                    target_path_length=target_path_length,
                    target_crow_length=target_crow_length,
                    internal_angle_mod=0.0,
                    N_turns=11,
                    mirrored=False,
                    resolution=200
                )
            else:
                wgs[idx].add_straight_segment(target_crow_length)
            wgs[idx].add_bend(-np.pi/2, param['min_radius'])
            wgs[idx].add_straight_segment((3-idx) * wg_sep - 50)

    x_extrema = (x_middle
                 - param['wg_sep']/2
                 - param['electrode_wg_sep'],
                 first_line_right_bound)
    wf_line_bounds = (x_extrema[0] - param['wf_leeways'][0],
                      initial_position_y,
                      x_extrema[1]
                      + param['wf_leeways'][0],
                      wgs[0].current_port.origin[1])
    line_bounds, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=wf_line_bounds,
        region_marker=region_marker,
        wf_maxlength=param['wf_maxlength'],
        wf_layer=param['wf_layer']
    )
    wg_wf_bounds = update_bounds(wg_wf_bounds, line_bounds)
    wf_line_bounds = [wf_line_bounds[0],
                      wf_line_bounds[3],
                      (wgs[0].current_port.origin[0]
                       + param['wf_leeways'][0]),
                      (wgs[0].get_shapely_outline().bounds[3]
                       + param['wf_leeways'][1])]
    top_line_bounds, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=wf_line_bounds,
        region_marker=region_marker,
        wf_maxlength=(wf_line_bounds[2] - wf_line_bounds[0])/3,
        wf_layer=param['wf_layer'],
        axis=0
    )
    wgs = wgs[::-1]
    wg_wf_bounds = update_bounds(wg_wf_bounds, top_line_bounds)
    initial_position_y = wgs[0].current_port.origin[1]
    # Second section 1
    for electrode_idx in range(2):
        sec1_electrodes = section_1_dcs_electrodes(
            cell=cell,
            wgs=wgs,
            electrode_wf_layers=param['right_electrode_wf_layers'][
                2*electrode_idx:2*(electrode_idx+2)
            ],
            ending_taper=(electrode_idx+1) % 2,
            param=param
        )
        right_side_electrodes[2*electrode_idx] = sec1_electrodes[0]
        right_side_electrodes[2*electrode_idx+1] = sec1_electrodes[1]
    # Second section 2
    # Expand the wgs first
    _, _ = expand_wgs_section_2(wgs, param)
    for electrode_idx in range(4, 6):
        sec2_electrodes = section_2_dcs_electrodes(
            cell=cell,
            wgs=wgs,
            electrode_wf_layer=(
                param['right_electrode_wf_layers'][electrode_idx]
            ),
            ending_taper=(electrode_idx+1) % 2,
            delay=delay,
            param=param
        )
        right_side_electrodes[electrode_idx] = sec2_electrodes
    
    # need to add the straight segments if there are no delay lines
    if not delay:
        wgs[0].add_straight_segment_until_y(wgs[1].current_port.origin[1])
        wgs[3].add_straight_segment_until_y(wgs[1].current_port.origin[1])
   
    # Add a second line of writefields following the waveguides
    second_line_x_diff = (wgs[-1].get_shapely_outline().bounds[2]
                          - second_x_middle)
    x_extrema = (second_x_middle
                 - second_line_x_diff,
                 second_x_middle
                 + second_line_x_diff)
    wf_line_bounds = (x_extrema[0] - param['wf_leeways'][0],
                      wgs[0].current_port.origin[1],
                      x_extrema[1]
                      + param['wf_leeways'][0],
                      top_line_bounds[1])
    line_bounds, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=wf_line_bounds,
        region_marker=region_marker,
        wf_maxlength=param['wf_maxlength'],
        wf_layer=param['wf_layer']
    )
    wg_wf_bounds = update_bounds(wg_wf_bounds, line_bounds)
    for idx, wg in enumerate(wgs):
        cell.add_to_layer(param['wg_layer'], wg)
    kink = [-1, -1]
    for idx in range(3, -1, -1):
        return_values = left_side_electrodes[idx]
        short_wg, signal_wg, long_wg = return_values[:3]
        direction, wf_line_bounds = return_values[3:]
        if direction == 1:
            kink = electrode_connector(
                cell=cell,
                left_wg=short_wg,
                signal_wg=signal_wg,
                right_wg=long_wg,
                kink=kink,
                wg_top=top_line_bounds[3],
                electrode_wf_layer=param['left_electrode_wf_layers'][idx],
                region_marker=left_electrode_region_markers[idx],
                connector_coordinates=left_connector_coordinates[idx],
                previous_line_bounds=wf_line_bounds,
                left_side=True,
                initial_point=initial_wf_point,
                param=param
            )
        else:
            kink = electrode_connector(
                cell=cell,
                left_wg=long_wg,
                signal_wg=signal_wg,
                right_wg=short_wg,
                kink=kink,
                wg_top=top_line_bounds[3],
                electrode_wf_layer=param['left_electrode_wf_layers'][idx],
                region_marker=left_electrode_region_markers[idx],
                connector_coordinates=left_connector_coordinates[idx],
                previous_line_bounds=wf_line_bounds,
                left_side=True,
                initial_point=initial_wf_point,
                param=param
            )
    kink = [-1, -1]
    for idx in range(6):
        return_values = right_side_electrodes[idx]
        short_wg, signal_wg, long_wg = return_values[:3]
        direction, wf_line_bounds = return_values[3:]
        if direction == 1:
            kink = electrode_connector(
                cell=cell,
                left_wg=long_wg,
                signal_wg=signal_wg,
                right_wg=short_wg,
                kink=kink,
                wg_top=top_line_bounds[3],
                electrode_wf_layer=param['right_electrode_wf_layers'][idx],
                region_marker=right_electrode_region_markers[idx],
                connector_coordinates=right_connector_coordinates[idx],
                previous_line_bounds=wf_line_bounds,
                left_side=False,
                initial_point=initial_wf_point,
                param=param
            )
        else:
            kink = electrode_connector(
                cell=cell,
                left_wg=short_wg,
                signal_wg=signal_wg,
                right_wg=long_wg,
                kink=kink,
                wg_top=top_line_bounds[3],
                electrode_wf_layer=param['right_electrode_wf_layers'][idx],
                region_marker=right_electrode_region_markers[idx],
                connector_coordinates=right_connector_coordinates[idx],
                previous_line_bounds=wf_line_bounds,
                left_side=False,
                initial_point=initial_wf_point,
                param=param
            )
    return wgs, wg_wf_bounds, region_marker


def wgs_to_fiber_array(
    cell,
    ports,
    coupler_positions,
    region_marker,
    is_incoupling,
    param,
    all_wf_layers=[],
    all_region_markers=[],
    alignment_test=False
):
    total_bounds = (np.inf, np.inf, -np.inf, -np.inf)
    first_line_bounds = list(total_bounds).copy()
    wgs = [Waveguide.make_at_port(prt) for prt in ports]
    first_bend = -np.pi/2
    if is_incoupling:
        first_bend = np.pi/2
    coupler_positions = coupler_positions[::-1]
    for idx, wg in enumerate(wgs):
        if is_incoupling:
            wg._current_port.angle = -np.pi/2
            wg.add_straight_segment(param['wg_sep']*(3-idx))
        else:
            wg.add_straight_segment(param['wg_sep']*(idx))
        # first bend then heading away from structure
        wg.add_bend(first_bend, param['min_radius'])
        if is_incoupling:
            wg.add_straight_segment_until_x(
                coupler_positions[-1][0] - param['min_radius']
            )
            if idx == 0:
                first_line_bounds = list(
                    update_bounds(
                        first_line_bounds,
                        wg.get_shapely_outline().bounds
                    )
                )
                # if alignment_test:
                first_line_bounds[2] -= 2*127
                first_line_bounds[0] -= 10
                first_line_bounds[1] -= 10
                if region_marker is None:
                    region_marker = bounds_to_polygon(first_line_bounds)
                _, region_marker = wf_line_from_bounds(
                    cell=cell,
                    bounds=first_line_bounds,
                    region_marker=region_marker,
                    wf_maxlength=1040,
                    wf_layer=param['wf_layer'],
                    axis=0,
                    direction=1
                )

            wg.add_straight_segment_until_x(
                coupler_positions[idx][0] - param['min_radius']
            )
            wg.add_bend(first_bend, param['min_radius'])
            wg.add_straight_segment_until_y(
                coupler_positions[idx][1]
            )
            gc = GratingCoupler.make_traditional_coupler_at_port(
                wg.current_port, **param['gc_params']
            )
            if idx == 0:
                second_line_bounds = [
                    first_line_bounds[2],
                    first_line_bounds[1],
                    (wg.current_port.origin[0]
                     + 127/2 - 10),
                    first_line_bounds[1] + 1040
                ]
            cell.add_to_layer(param['wg_layer'], gc, wg)
            add_markers_to_top(
                cell=cell,
                wf_bounds=second_line_bounds,
                device_top_bound=(
                    wgs[0].current_port.origin[1]+70
                ),
                marker_dims=20,
                wg_layer=param['wg_layer'],
                marker_layer_1=3,
                marker_layer_2=4,
                marker_protection_layer=15
            )

        else:
            wg.add_straight_segment_until_x(
                coupler_positions[0][0] + param['min_radius']
            )
            if idx == 3:
                first_line_bounds = list(update_bounds(
                    first_line_bounds,
                    wg.get_shapely_outline().bounds
                ))
                first_line_bounds[1] -= 10
                first_line_bounds[2] += 10
                # if alignment_test:
                first_line_bounds[0] += 1.5*127
                if region_marker is None:
                    region_marker = bounds_to_polygon(first_line_bounds)
                _, region_marker = wf_line_from_bounds(
                    cell=cell,
                    bounds=first_line_bounds,
                    region_marker=region_marker,
                    wf_maxlength=1040,
                    wf_layer=param['wf_layer'],
                    axis=0,
                    direction=-1
                )

            wg.add_straight_segment_until_x(
                coupler_positions[idx][0] + param['min_radius']
            )
            wg.add_bend(first_bend, param['min_radius'])
            wg.add_straight_segment_until_y(
                coupler_positions[idx][1]
            )
            gc = GratingCoupler.make_traditional_coupler_at_port(
                wg.current_port, **param['gc_params']
            )
            cell.add_to_layer(param['wg_layer'], gc, wg)
            if idx == 3:
                second_line_bounds = [
                    (wg.current_port.origin[0]
                     - 127/2+10),
                    first_line_bounds[1],
                    first_line_bounds[0],
                    first_line_bounds[1] + 1040
                ]
    if alignment_test:
        test_coupler_positions = [
            [
                coupler_positions[-1][0] - 127,
                coupler_positions[-1][1]
            ],
            [
                coupler_positions[-1][0] - 2*127,  # +8*127
                coupler_positions[-1][1]
            ]
        ]
        test_port = Port(test_coupler_positions[0],
                         np.pi/2.0,
                         wgs[0].width)
        test_gc_0 = GratingCoupler.make_traditional_coupler_at_port(
            test_port, **param['gc_params']
        )
        test_wg = Waveguide.make_at_port(test_gc_0.port)
        test_wg._current_port.angle = -np.pi/2
        test_wg.add_bend(-np.pi/2, param['min_radius'])
        test_wg.add_straight_segment_until_x(
            test_coupler_positions[1][0] + param['min_radius']
        )
        test_wg.add_bend(-np.pi/2, param['min_radius'])
        test_gc_1 = GratingCoupler.make_traditional_coupler_at_port(
            test_wg.current_port, **param['gc_params']
        )
        cell.add_to_layer(param['wg_layer'], test_wg, test_gc_0, test_gc_1)

    _, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=second_line_bounds,
        region_marker=region_marker,
        wf_maxlength=1040,
        wf_layer=param['wf_layer'],
        axis=1,
        direction=1
    )
    if is_incoupling:
        for idx in range(len(all_wf_layers)):
            all_region_markers[idx] = bounds_to_polygon(second_line_bounds)
            _, all_region_markers[idx] = wf_line_from_bounds(
                cell=cell,
                bounds=second_line_bounds,
                region_marker=all_region_markers[idx],
                wf_maxlength=1040,
                wf_layer=all_wf_layers[idx],
                axis=1,
                direction=1
            )

    total_bounds = update_bounds(first_line_bounds, second_line_bounds)
    return total_bounds, region_marker


def interferometer_and_fiber_array(
    cell,
    inports,
    gc_positions,
    electrode_contact_pad_coordinates,
    second_x_middle,
    delay_lines=True,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    **kwargs
):
    DEFAULT_PARAMETERS = {
        'sm_wg_width': 0.5,
        'mm_wg_width': 0.5,
        'mm_taper_length': 0,
        'min_radius': 50,
        'mzi_sep_leeway': 50,
        'wg_sep': 25,
        'electrode_wg_sep': 100,
        'last_layer': False,
        'delay_line_length': 100,
        # For directional couplers and bends
        'sine_s_x': 60,
        'coupler_sep': 0.4,
        'coupler_length': 18,
        # For electrodes
        'electrode_length': 1250,
        # For electrode function
        'electrode_width': 25,
        'electrode_sep': 1.1,
        'crossing_width': 10,
        'electrode_taper_leeway': 5,
        'electrode_taper_length': 30,
        'electrode_sep_y': 15,
        'electrode_pitch': 40,
        'connector_start_y': 8000,
        'connector_probe_dims': (80, 600),
        'connector_probe_pitch': 150,
        'contact_pad_dims': (250, 600),
        'contact_pad_sep': 150,
        'electrode_pad_seps': (30, 30),
        'pad_sep': 300,
        # For grating couplers
        'alignment_test': True,
        # For layers and writefields
        'wg_layer': 9,
        'wf_layer': 100,
        'electrode_layer': 10,
        'right_electrode_wf_layers': np.arange(110, 140+1, 10),
        'left_electrode_wf_layers': np.arange(150, 200+1, 10),
        # For writefields
        'wf_maxlength': 1040,
        'wf_leeways': (10, 10),
        'wf_x_sep': 10,
        'wf_electrode_leeways': (10, 10),
        'wg_region_layer': 101,
        # For markers
        'marker_dims': 20,
        'mk_layer_1': 3,
        'mk_layer_2': 4,
        'mk_layer_3': 15,
        # For GC
        'gc_params': {
            'width': 0.5,
            'full_opening_angle': np.deg2rad(180),
            'grating_period': 0.46,
            'grating_ff': 0.2,
            'n_gratings': 10,
            'ap_max_ff': 0.8,
            'n_ap_gratings': 65,
            'taper_length': 9.64
        }
    }

    # From other function
    param = DEFAULT_PARAMETERS.copy()
    fix_dict(param, kwargs)
    wgs = [Waveguide.make_at_port(inport) for inport in inports]
    region_marker = None
    electrode_region_markers = (
        [None] * (len(param['right_electrode_wf_layers'])
                  + len(param['left_electrode_wf_layers']))
    )
    electrode_wfs = np.concatenate((param['right_electrode_wf_layers'],
                                    param['left_electrode_wf_layers']))
    incoupling_bounds, region_marker = wgs_to_fiber_array(
        cell=cell,
        ports=inports,
        coupler_positions=gc_positions[:4],
        region_marker=region_marker,
        is_incoupling=True,
        all_wf_layers=electrode_wfs,
        all_region_markers=electrode_region_markers,
        param=param
    )
    right_electrode_region_markers = (
        electrode_region_markers[:len(param['right_electrode_wf_layers'])]
    )
    left_electrode_region_markers = (
        electrode_region_markers[len(param['right_electrode_wf_layers']):]
    )

    initial_wf_point = ((incoupling_bounds[0] + incoupling_bounds[2])/2,
                        incoupling_bounds[3])

    if delay_lines:
        wgs, total_bounds, region_marker = build_interferometer(
            cell,
            wgs,
            right_electrode_region_markers,
            left_electrode_region_markers,
            region_marker,
            second_x_middle,
            initial_wf_point,
            delay=True,
            param=param,
            total_bounds=total_bounds,
        )
    else:
        wgs, total_bounds, region_marker = build_interferometer(
            cell,
            wgs,
            right_electrode_region_markers,
            left_electrode_region_markers,
            region_marker,
            second_x_middle,
            initial_wf_point,
            delay=False,
            param=param,
            total_bounds=total_bounds,
        )
    outports = [wg.current_port for wg in wgs]

    total_bounds = update_bounds(total_bounds, incoupling_bounds)

    outcoupling_bounds, region_marker = wgs_to_fiber_array(
        cell=cell,
        ports=outports,
        coupler_positions=gc_positions[4:],
        region_marker=region_marker,
        param=param,
        is_incoupling=False
    )

    cell.add_to_layer(param['wg_region_layer'], region_marker)
    total_bounds = update_bounds(total_bounds, outcoupling_bounds)

    return total_bounds
