import numpy as np
from shapely.geometry import Polygon

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


def wf_line_from_bounds(cell, bounds, wf_maxlength, wf_layer, wf_additions = {}, axis=1):
    # assumes x direction is fine and that leeways have been included in the bounds
    current_position = bounds[axis]
    current_limit = min(current_position + wf_maxlength, bounds[axis+2])
    wf_bounds = list(bounds).copy()
    iteration = 0
    while current_limit != bounds[axis+2]:
        wf_bounds[axis] = current_position
        wf_bounds[axis+2] = current_limit
        wf_addition = (0, 0, 0, 0)
        if iteration in wf_additions.keys():
            wf_addition = wf_additions[iteration]
            bounds = np.add(bounds, wf_addition)
        single_wf_from_bounds(cell, np.add(wf_bounds, wf_addition), wf_layer)
        current_position=current_limit
        current_limit = min(current_position + wf_maxlength, bounds[axis+2])
        iteration+=1
    wf_addition = (0, 0, 0, 0)
    if -1 in wf_additions.keys():
        wf_addition = wf_additions[-1]
        bounds = np.add(bounds, wf_addition)
    wf_bounds[axis] = current_position
    wf_bounds[axis+2] = current_limit
    single_wf_from_bounds(cell, np.add(wf_bounds, wf_addition), wf_layer)
    return bounds
