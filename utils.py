import numpy as np
from shapely.geometry import Polygon
from gdshelpers.geometry.shapely_adapter import geometric_union
from gdshelpers.parts.marker import DLWMarker, SquareMarker, CrossMarker

def update_bounds(bounds1, bounds2):
    # (minx, miny, maxx, maxy)
    return (min(bounds1[0], bounds2[0]),
            min(bounds1[1], bounds2[1]),
            max(bounds1[2], bounds2[2]),
            max(bounds1[3], bounds2[3]))


def fix_dict(parameters, kwargs):
    if 'kwarg_verbose' in kwargs:
        parameters['kwarg_verbose'] = kwargs['kwarg_verbose']
    else:
        parameters['kwarg_verbose'] = False
    for key, value in kwargs.items():
        if key in parameters:
            parameters[key] = value
        else:
            if parameters['kwarg_verbose']:
                print(key, 'is not a parameter for this function')


def single_wf_from_bounds(cell, wf_bounds, wf_layer):
    wf_corners = [(  # top_left corner (minx, maxy)
                          wf_bounds[0],
                          wf_bounds[3]
                      ),
                      (  # top_right corner (maxx, maxy)
                          wf_bounds[2],
                          wf_bounds[3]
                      ),
                      (  # bottom_right corner (maxx, miny)
                          wf_bounds[2],
                          wf_bounds[1]
                      ),
                      (  # bottom_left_corner (minx, miny)
                          wf_bounds[0],
                          wf_bounds[1]
                      )]
    wf_polygon = Polygon(wf_corners)
    cell.add_to_layer(wf_layer, wf_polygon)


def wf_line_from_bounds(cell,
                        bounds,
                        wf_maxlength,
                        wf_layer,
                        wf_additions={},
                        axis=1,
                        region_marker=None,
                        direction=1):
    # assumes x direction is fine and that leeways have
    # been included in the bounds
    dir_idx = 1
    if direction == -1:
        dir_idx = 0
    current_position = bounds[axis+2*(1-dir_idx)]
    if direction == 1:
        current_limit = min(current_position + wf_maxlength, bounds[axis+2])
    else:
        current_limit = max(current_position - wf_maxlength, bounds[axis])
    wf_bounds = list(bounds).copy()
    iteration = 0
    while current_limit != bounds[axis+2*dir_idx]:
        wf_bounds[axis+2*(1-dir_idx)] = current_position
        wf_bounds[axis+2*dir_idx] = current_limit
        wf_addition = (0, 0, 0, 0)
        if iteration in wf_additions.keys():
            wf_addition = wf_additions[iteration]
            bounds = np.add(bounds, wf_addition)
        single_wf_from_bounds(cell, np.add(wf_bounds, wf_addition), wf_layer)
        current_position = current_limit
        if direction == 1:
            current_limit = min(current_position + wf_maxlength,
                                bounds[axis+2])
        else:
            current_limit = max(current_position - wf_maxlength, bounds[axis])
        iteration += 1
    wf_addition = (0, 0, 0, 0)
    if -1 in wf_additions.keys():
        wf_addition = wf_additions[-1]
        bounds = np.add(bounds, wf_addition)
    wf_bounds[axis+2*(1-dir_idx)] = current_position
    wf_bounds[axis+2*dir_idx] = current_limit
    single_wf_from_bounds(cell, np.add(wf_bounds, wf_addition), wf_layer)
    if region_marker is not None:
        region_marker = geometric_union([
            region_marker,
            bounds_to_polygon(bounds)
        ])
    return bounds, region_marker


def bounds_to_polygon(bounds):
    # (minx, maxy), (maxx, maxy), (maxx, miny), (minx, miny)
    corners = [(bounds[0], bounds[3]), (bounds[2], bounds[3]), 
               (bounds[2], bounds[1]), (bounds[0], bounds[1])]
    return Polygon(corners)


def add_markers_to_top(cell,
                       wf_bounds,
                       device_top_bound,
                       marker_dims,
                       wg_layer,
                       marker_layer_1,
                       marker_layer_2,
                       marker_protection_layer):
    marker_positions_1 = [(wf_bounds[0] + 1.5*marker_dims+20,
                           wf_bounds[3] - 1.5*marker_dims-80),
                          (wf_bounds[2] - 1.5*marker_dims,
                           wf_bounds[3] - 1.5*marker_dims-80),
                          (wf_bounds[2] - 1.5*marker_dims,
                           device_top_bound + 80
                           + 3*marker_dims)]
    marker_positions_2 = [(marker_positions_1[0][0] + 120,
                           marker_positions_1[0][1]),
                          (marker_positions_1[1][0] - 120,
                           marker_positions_1[1][1]),
                          (marker_positions_1[2][0] - 120,
                           marker_positions_1[2][1])]
    markers_1 = [SquareMarker.make_marker(position, marker_dims)
                 for position in marker_positions_1]
    markers_1_wg = [SquareMarker.make_marker(position, 1.5*marker_dims)
                    for position in marker_positions_1]
    markers_1_protection = [SquareMarker.make_marker(position, 2*marker_dims)
                            for position in marker_positions_1]
    markers_2 = [SquareMarker.make_marker(position, marker_dims)
                 for position in marker_positions_2]
    markers_2_wg = [SquareMarker.make_marker(position, 1.5*marker_dims)
                    for position in marker_positions_2]
    markers_2_protection = [SquareMarker.make_marker(position, 2*marker_dims)
                            for position in marker_positions_2]
    for idx in range(len(markers_1)):
        cell.add_to_layer(marker_layer_1, markers_1[idx])
        cell.add_to_layer(wg_layer, markers_1_wg[idx])
        cell.add_to_layer(marker_protection_layer,
                          markers_1_protection[idx])
        cell.add_to_layer(marker_layer_2, markers_2[idx])
        cell.add_to_layer(wg_layer, markers_2_wg[idx])
        cell.add_to_layer(marker_protection_layer,
                          markers_2_protection[idx])


def connect_writefields(cell,
                        initial_point,
                        end_point,
                        wf_layer,
                        region_marker,
                        parameters):
    if end_point[1] < initial_point[1]:
        new_end_point = [end_point[0],
                         initial_point[1] + (initial_point[1]-end_point[1])]
        end_point = new_end_point

    vertical_line_bounds = [
                            initial_point[0]-parameters['wf_leeways'][0],
                            initial_point[1],
                            initial_point[0]+parameters['wf_leeways'][0],
                            end_point[1]-parameters['wf_leeways'][0]
                           ]
    _, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=vertical_line_bounds,
        region_marker=region_marker,
        wf_maxlength=1040,
        wf_layer=wf_layer,
        axis=1,
        direction=1
    )
    if end_point[0] > initial_point[0]:
        horizontal_line_bounds = [
                                  initial_point[0]-parameters['wf_leeways'][0],
                                  end_point[1]-parameters['wf_leeways'][0],
                                  end_point[0],
                                  end_point[1]+parameters['wf_leeways'][0]
                                 ]
    else:
        horizontal_line_bounds = [
                                  end_point[0],
                                  end_point[1]-parameters['wf_leeways'][0],
                                  initial_point[0]+parameters['wf_leeways'][0],
                                  end_point[1]+parameters['wf_leeways'][0]
                                 ]
    _, region_marker = wf_line_from_bounds(
        cell=cell,
        bounds=horizontal_line_bounds,
        region_marker=region_marker,
        wf_maxlength=1040,
        wf_layer=wf_layer,
        axis=1,
        direction=1
    )
    return region_marker
