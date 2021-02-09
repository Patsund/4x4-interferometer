import numpy as np
from gdshelpers.geometry.chip import Cell
from gdshelpers.parts.text import Text
from gdshelpers.parts.marker import DLWMarker, SquareMarker, CrossMarker
from gdshelpers.parts.waveguide import Waveguide
from gdshelpers.parts.port import Port
from gdshelpers.geometry.shapely_adapter import geometric_union
from gdshelpers.parts.coupler import GratingCoupler
from gdshelpers.parts.spiral import Spiral
from shapely.geometry import Polygon
from utils import *
from parts import build_interferometer, wgs_to_fiber_array
from parts import add_markers_to_bottom_of_cell

cell = Cell('MZI_active')

wg_width = 0.5
wafer_min_max = (0, 7200)

electrode_wg_sep = 62.5

inport_0 = Port((-electrode_wg_sep/2 - 25, 0), np.pi/2, wg_width)
inport_1 = Port((-electrode_wg_sep/2, 0), np.pi/2, wg_width)
inport_2 = Port((electrode_wg_sep/2, 0), np.pi/2, wg_width)
inport_3 = Port((electrode_wg_sep/2 + 25, 0), np.pi/2, wg_width)
device_inports = [inport_0, inport_1, inport_2, inport_3]
wgs_1 = [Waveguide.make_at_port(inport) for inport in device_inports]

wgs, _ = build_interferometer(
    cell=cell,
    wgs=wgs_1,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=9,
    wf_layer=100,
    electrode_layer=10,
    electrode_wf_layer=102,
    electrode_contact_pad_coordinates=[],
    second_x_middle=2 * 127,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width=1,
    mm_taper_length=30,
    electrode_wg_sep=electrode_wg_sep
)

# fiber_array
gc_pitch = 127
gc_leeway = 100
wg_sep = 25
min_radius = 50
starting_point = -gc_pitch * 2.5
ports = [wg.current_port for wg in wgs]
gc_y = (device_inports[0].origin[1]
        - min_radius
        - 3*wg_sep
        - gc_leeway)
grating_coupler_positions = [(starting_point+idx*gc_pitch, gc_y)
                             for idx in range(8)]


incoupling_bounds = wgs_to_fiber_array(
    cell=cell,
    ports=device_inports,
    coupler_positions=grating_coupler_positions[:4],
    min_radius=min_radius,
    gc_leeway=gc_leeway,
    wg_sep=wg_sep,
    wg_layer=9,
    wf_layer=100,
    is_incoupling=True
)

incoupling_wf_bounds = [incoupling_bounds[0]-15,
                        incoupling_bounds[3]-690,
                        incoupling_bounds[2]+gc_pitch/2-15,
                        incoupling_bounds[3]]

_ = wf_line_from_bounds(
            cell=cell,
            bounds=incoupling_wf_bounds,
            wf_maxlength=1040,
            wf_layer=100,
            axis=1
)

add_markers_to_bottom_of_cell(cell=cell,
                              wf_bounds=incoupling_wf_bounds,
                              device_bottom_bound=incoupling_bounds[1],
                              marker_dims=20,
                              marker_layer_1=3,
                              marker_layer_2=4,
                              marker_protection_layer=15)

outcoupling_bounds = wgs_to_fiber_array(
    cell=cell,
    ports=ports,
    coupler_positions=grating_coupler_positions[4:],
    min_radius=min_radius,
    gc_leeway=gc_leeway,
    wg_sep=wg_sep,
    wg_layer=9,
    wf_layer=100,
    is_incoupling=False
)

outcoupling_wf_bounds = [outcoupling_bounds[0]-gc_pitch/2+15,
                         outcoupling_bounds[1]-15,
                         outcoupling_bounds[2]+15,
                         outcoupling_bounds[3]]

_ = wf_line_from_bounds(
            cell=cell,
            bounds=outcoupling_wf_bounds,
            wf_maxlength=1040,
            wf_layer=100,
            axis=1
)

cell.save('4x4_secondgen.gds')
cell.show()
