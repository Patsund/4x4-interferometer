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
from parts import interferometer_and_fiber_array

cell = Cell('4x4 interferometer')

wg_width = 0.5
electrode_wg_sep = 62.5
second_x_middle = 8 * 2 * 127

inport_0 = Port((-electrode_wg_sep/2 - 25, 0), np.pi/2, wg_width)
inport_1 = Port((-electrode_wg_sep/2, 0), np.pi/2, wg_width)
inport_2 = Port((electrode_wg_sep/2, 0), np.pi/2, wg_width)
inport_3 = Port((electrode_wg_sep/2 + 25, 0), np.pi/2, wg_width)
device_inports = [inport_0, inport_1, inport_2, inport_3]

# fiber_array
gc_pitch = 127
wg_sep = 25
min_radius = 50
starting_point = second_x_middle/2 - gc_pitch * 3.5
gc_y = 100
grating_coupler_positions = [(starting_point+idx*gc_pitch, gc_y)
                             for idx in range(8)]


_ = interferometer_and_fiber_array(
    cell=cell,
    inports=device_inports,
    gc_positions=grating_coupler_positions,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=9,
    wf_layer=100,
    electrode_layer=10,
    electrode_wf_layer=102,
    electrode_contact_pad_coordinates=[],
    second_x_middle=second_x_middle,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width=1,
    mm_taper_length=30,
    contact_pad_dims=[250, 250],
    contact_pad_sep=160,
    connector_start_y=7120,
    connector_probe_pitch=150,
    wg_sep=wg_sep,
    electrode_wg_sep=electrode_wg_sep
)

cell.save('4x4_3rdgen_150_1.gds')
