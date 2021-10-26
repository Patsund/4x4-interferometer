import numpy as np
from utils import *
from gdshelpers.parts.waveguide import Waveguide
from gdshelpers.parts.port import Port
from gdshelpers.parts.coupler import GratingCoupler


def expand_wgs_section_2(wgs, param):
    # assume the first (final) inport is always the leftmost (rightmost)
    global_x_middle = (wgs[0].current_port.origin[0]
                       + wgs[-1].current_port.origin[0])/2
    desired_x_position = [global_x_middle + (idx - 1.5) * param['wg_sep']
                          for idx in range(len(wgs))]
    desired_x_position[0] = (desired_x_position[1]
                             - param['electrode_wg_sep'])
    desired_x_position[3] = (desired_x_position[2]
                             + param['electrode_wg_sep'])
    # get the wgs to where we want them to be
    for idx, wg in enumerate(wgs):
        if wg.current_port.angle == np.pi/2:
            s_y_length = desired_x_position[idx] - wg.current_port.origin[0]
        else:
            s_y_length = wg.current_port.origin[0] - desired_x_position[idx]
        s_x_length = param['sine_s_x']
        wg.add_parameterized_path(path=lambda t: (t * s_x_length, .5
                                                  * (np.cos(np.pi * t) - 1)
                                                  * s_y_length),
                                  path_derivative=lambda t: (
                                      s_x_length, -np.pi * .5
                                      * np.sin(np.pi * t) * s_y_length
                                  ))
    return global_x_middle, desired_x_position


def just_dc_wg(
    wgs,
    param,
    stagger_separation=0,
    ending_taper=True
):
    for idx, wg in enumerate(wgs):
        increase_length = False
        if (
            wg.current_port.width < param['mm_wg_width']
            and param['mm_taper_length'] > 0
        ):
            wg.add_straight_segment(param['mm_taper_length'],
                                    param['mm_wg_width'])
            increase_length = True

        # taper for mm assuming the wg is mm during electrodes
        if param['mm_taper_length'] > 0:
            wg.add_straight_segment(param['mm_taper_length'],
                                    param['sm_wg_width'])

        # directional coupler
        x_length = param['sine_s_x']
        y_length = (param['wg_sep'] / 2.0
                    - (param['coupler_sep'] + param['sm_wg_width']) / 2.0
                    - idx * (param['wg_sep'] - param['sm_wg_width']
                             - param['coupler_sep'])
                    )
        if wg.current_port.angle == -np.pi/2:
            y_length *= -1
        wg.add_parameterized_path(path=lambda t: (t * x_length, .5
                                                  * (np.cos(np.pi * t) - 1)
                                                  * y_length),
                                  path_derivative=lambda t: (
                                         x_length, -np.pi * .5
                                         * np.sin(np.pi * t) * y_length
                                  ))
        wg.add_straight_segment(param['coupler_length'])
        y_length = -(param['wg_sep'] / 2.0
                     - (param['coupler_sep'] + param['sm_wg_width']) / 2.0
                     - idx * (param['wg_sep']
                     - param['sm_wg_width']
                     - param['coupler_sep']))
        if wg.current_port.angle == -np.pi/2:
            y_length *= -1
        wg.add_parameterized_path(path=lambda t: (t * x_length, .5
                                                  * (np.cos(np.pi * t) - 1)
                                                  * y_length),
                                  path_derivative=lambda t: (
                                      x_length, -np.pi * .5
                                      * np.sin(np.pi * t) * y_length
                                  ))
        if param['mm_taper_length'] > 0 and ending_taper:
            wg.add_straight_segment(param['mm_taper_length'],
                                    param['mm_wg_width'])

def phase_shifter_and_dc_wg(
    wgs,
    param,
    stagger_separation=0,
    ending_taper=True
):
    for idx, wg in enumerate(wgs):
        increase_length = False
        # straight segment for the electrode
        if (
            wg.current_port.width < param['mm_wg_width']
            and param['mm_taper_length'] > 0
        ):
            wg.add_straight_segment(param['mm_taper_length'],
                                    param['mm_wg_width'])
            increase_length = True
        wg.add_straight_segment(param['electrode_length']+stagger_separation)

        # taper for mm assuming the wg is mm during electrodes
        if param['mm_taper_length'] > 0:
            wg.add_straight_segment(param['mm_taper_length'],
                                    param['sm_wg_width'])

        # directional coupler
        x_length = param['sine_s_x']
        y_length = (param['wg_sep'] / 2.0
                    - (param['coupler_sep'] + param['sm_wg_width']) / 2.0
                    - idx * (param['wg_sep'] - param['sm_wg_width']
                             - param['coupler_sep'])
                    )
        if wg.current_port.angle == -np.pi/2:
            y_length *= -1
        wg.add_parameterized_path(path=lambda t: (t * x_length, .5
                                                  * (np.cos(np.pi * t) - 1)
                                                  * y_length),
                                  path_derivative=lambda t: (
                                         x_length, -np.pi * .5
                                         * np.sin(np.pi * t) * y_length
                                  ))
        wg.add_straight_segment(param['coupler_length'])
        y_length = -(param['wg_sep'] / 2.0
                     - (param['coupler_sep'] + param['sm_wg_width']) / 2.0
                     - idx * (param['wg_sep']
                     - param['sm_wg_width']
                     - param['coupler_sep']))
        if wg.current_port.angle == -np.pi/2:
            y_length *= -1
        wg.add_parameterized_path(path=lambda t: (t * x_length, .5
                                                  * (np.cos(np.pi * t) - 1)
                                                  * y_length),
                                  path_derivative=lambda t: (
                                      x_length, -np.pi * .5
                                      * np.sin(np.pi * t) * y_length
                                  ))
        if param['mm_taper_length'] > 0 and ending_taper:
            wg.add_straight_segment(param['mm_taper_length'],
                                    param['mm_wg_width'])


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
            # adiabatic taper
            wg.add_straight_segment(40, final_width=0.75)
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
            # adiabatic taper
            wg.add_straight_segment(40, final_width=0.75)
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
