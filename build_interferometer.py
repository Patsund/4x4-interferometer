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
from parts import build_interferometer
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
    wg_layer=7,
    wf_layer=107,
    electrode_layer=15,
    electrode_wf_layer=115,
    electrode_contact_pad_coordinates=[],
    second_x_middle=2 * 127,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width=1,
    mm_taper_length=30,
    electrode_wg_sep=electrode_wg_sep
)

cell.save('4x4_secondgen_test.gds')
cell.show()