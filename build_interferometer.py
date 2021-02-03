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
from parts import build_section

cell = Cell('MZI_active')

wg_width = 0.5
wafer_min_max = (0, 7200)

inport_0 = Port((-25, 0), np.pi/2, wg_width)
inport_1 = Port((0, 0), np.pi/2, wg_width)
inport_2 = Port((25, 0), np.pi/2, wg_width)
inport_3 = Port((2*25, 0), np.pi/2, wg_width)
device_inports = [inport_0, inport_1, inport_2, inport_3]

initial_positions = (-25, 0)
outports_1, initial_positions_1 = build_section(
    cell=cell,
    section=1,
    inports=device_inports,
    initial_positions = initial_positions,
    wafer_min_max=wafer_min_max,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=7,
    wf_layer=107,
    electrode_layer=15,
    electrode_wf_layer=115,
    electrode_sep=1.5,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width = 1,
    mm_taper_length = 30
)

outports_2, initial_positions_2 = build_section(
    cell=cell,
    section=2,
    inports=outports_1,
    initial_positions = initial_positions_1,
    wafer_min_max=wafer_min_max,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=7,
    wf_layer=107,
    electrode_layer=15,
    electrode_wf_layer=115,
    electrode_sep=1.5,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width = 1,
    mm_taper_length = 30
)

outports_3, initial_positions_3 = build_section(
    cell=cell,
    section=1,
    inports=outports_2,
    initial_positions = initial_positions_2,
    wafer_min_max=wafer_min_max,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=7,
    wf_layer=107,
    electrode_layer=15,
    electrode_wf_layer=115,
    electrode_sep=1.5,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width = 1,
    mm_taper_length = 30
)

outports_4, initial_positions_4 = build_section(
    cell=cell,
    section=2,
    inports=outports_3,
    initial_positions = initial_positions_3,
    wafer_min_max=wafer_min_max,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=7,
    wf_layer=107,
    electrode_layer=15,
    electrode_wf_layer=115,
    electrode_sep=1.5,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    last_layer = True,
    mm_wg_width = 1,
    mm_taper_length = 30
)

fiber_array_incouplers = [Waveguide.make_at_port(inport) for inport in device_inports]
fiber_array_outcouplers = [Waveguide.make_at_port(outport) for outport in outports_4]

wg_sep = 25
fiber_coupling_leeways = (100, 500)
fiber_coupling_y_position = -450
coupler_separation = 127.
min_radius = 50

for idx, wg in enumerate(fiber_array_incouplers):
    wg._current_port.angle = wg.current_port.angle + np.pi
    wg.add_bend(-np.pi/2, min_radius + wg_sep * idx)
    wg.add_straight_segment(300)
    wg.add_bend(np.pi/2, min_radius + wg_sep * (3-idx))
    # y_segment
    goal_y_position = fiber_coupling_y_position
    current_y_position = wg.current_port.origin[1]
    wg.add_straight_segment(current_y_position - goal_y_position - wg_sep*8)
    wg.add_bend(np.pi/2, min_radius + wg_sep*(3-idx))
    goal_x_position = ((device_inports[0].origin[0] + outports_4[-1].origin[0])/2 
                       - coupler_separation/2
                       - (3-idx) * coupler_separation - min_radius)
    current_x_position = wg.current_port.origin[0]
    wg.add_straight_segment(abs(goal_x_position-current_x_position))
    wg.add_bend(-np.pi/2, min_radius)
    current_y_position = wg.current_port.origin[1]
    wg.add_straight_segment(current_y_position-goal_y_position)
    grating_coupler = GratingCoupler.make_traditional_coupler_at_port(
        wg.current_port, **std_coupler_params
    )
    cell.add_to_layer(7, wg, grating_coupler)
    
for idx, wg in enumerate(fiber_array_outcouplers):
    wg.add_bend(np.pi/2, min_radius + wg_sep * (3-idx))
    wg.add_straight_segment(300)
    wg.add_bend(-np.pi/2, min_radius + wg_sep * idx)
    # y_segment
    goal_y_position = fiber_coupling_y_position
    current_y_position = wg.current_port.origin[1]
    wg.add_straight_segment(current_y_position - goal_y_position - wg_sep*8)
    wg.add_bend(-np.pi/2, min_radius + wg_sep*idx)
    goal_x_position = ((device_inports[0].origin[0] + outports_4[-1].origin[0])/2 
                       + coupler_separation/2
                       + (idx) * coupler_separation + min_radius)
    current_x_position = wg.current_port.origin[0]
    wg.add_straight_segment(current_x_position - goal_x_position)
    wg.add_bend(np.pi/2, min_radius)
    current_y_position = wg.current_port.origin[1]
    wg.add_straight_segment(current_y_position-goal_y_position)
    grating_coupler = GratingCoupler.make_traditional_coupler_at_port(
        wg.current_port, **std_coupler_params
    )
    cell.add_to_layer(7, wg, grating_coupler)

cell.save('4x4_secondgen_test.gds')
cell.show()
