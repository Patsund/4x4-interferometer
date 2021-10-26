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
from gdshelpers.helpers import id_to_alphanumeric
from utils import *
from parts import interferometer_and_fiber_array, build_global_markers
from TestStructures.tests_DC_Grating_Delays import RectangularSpiral, \
    DirectionalCouplersTest, Efficiency_Grating, Ring_Test, DirectionalCouplersTest_standard
from TestStructures.tests_Modulators import MZI_active
from TestStructures.test_TimeBinBS import TimeBin_BS
from TestStructures.tests_ModulatorsWithPhase import MZI_active_with_phase
from TestStructures.test_Demux import Demux_active

from tech.LiNb01 import *

wg_Expwidth = 0.5
std_coupler_params['grating_period'] = 0.49


electrode_wg_sep = 62.5
second_x_middle = 8 * 2 * 127

inport_0 = Port((-electrode_wg_sep / 2 - 25, 0), np.pi / 2, wg_width)
inport_1 = Port((-electrode_wg_sep / 2, 0), np.pi / 2, wg_width)
inport_2 = Port((electrode_wg_sep / 2, 0), np.pi / 2, wg_width)
inport_3 = Port((electrode_wg_sep / 2 + 25, 0), np.pi / 2, wg_width)
device_inports = [inport_0, inport_1, inport_2, inport_3]

start_x = inport_0.origin[0] + 3500  # should be + 3500 for 150 probe pitch
device_inports_2 = [None, None, None, None]
device_inports_2[:2] = [
    Port((start_x + idx * 25, 0), np.pi / 2, wg_width)
    for idx in range(2)
]
start_x += 25 + electrode_wg_sep
first_x_middle_2 = start_x - electrode_wg_sep / 2
second_x_middle_2 = first_x_middle_2 + second_x_middle
device_inports_2[2:] = [
    Port((start_x + idx * 25, 0), np.pi / 2, wg_width)
    for idx in range(2)
]

# fiber_array
gc_pitch = 127
wg_sep = 25
min_radius = 70
gc_starting_point = second_x_middle / 2 - gc_pitch * 3.5
gc_starting_point_2 = (
        0.5 * (first_x_middle_2 + second_x_middle_2)
        - gc_pitch * 3.5
)
gc_y = 100
grating_coupler_positions = [(gc_starting_point + idx * gc_pitch, gc_y)
                             for idx in range(8)]
grating_coupler_positions_2 = [(gc_starting_point_2 + idx * gc_pitch, gc_y)
                               for idx in range(8)]

# Adding Left 4x4 interferometer
left4x4_cell = Cell('4x4 interferometer_Left')

left_electrode_wf_layers = np.arange(110, 140+1, 10)
right_electrode_wf_layers = np.arange(150, 200+1, 10)


_, left_init_wf_point, left_in_bounds = interferometer_and_fiber_array(
    cell=left4x4_cell,
    inports=device_inports,
    gc_positions=grating_coupler_positions,
    electrode_length=1250,
    coupler_sep=0.5,
    coupler_length=30,
    sm_wg_width=wg_width,
    delay_lines=False,
    wg_layer=9,
    wf_layer=100,
    electrode_layer=10,
    right_electrode_wf_layers=right_electrode_wf_layers,
    left_electrode_wf_layers=left_electrode_wf_layers,
    electrode_contact_pad_coordinates=[],
    second_x_middle=second_x_middle,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width=wg_Expwidth,
    mm_taper_length=30,
    contact_pad_dims=[250, 250],
    contact_pad_sep=160,
    connector_start_y=7180,
    connector_probe_pitch=175,
    wg_sep=wg_sep,
    electrode_wg_sep=electrode_wg_sep,
    gc_params={
            'width': 0.5,
            'full_opening_angle': np.deg2rad(180),
            'grating_period': 0.49,
            'grating_ff': 0.2,
            'n_gratings': 10,
            'ap_max_ff': 0.8,
            'n_ap_gratings': 65,
            'taper_length': 9.64
    }
)
left4x4_cell.add_to_layer(comment_layer, Text((560, 540), 10,
                    'wg_width={:.2f}um\nexp_wg_width={:.2f}um\ngrat_period={:.2f}um'.format(wg_width,wg_Expwidth, 0.49),
                    alignment='center-top'))

# Adding Right 4x4 interferometer
right4x4_cell = Cell('4x4 interferometer_Right')
_, right_init_wf_point, right_in_bounds = interferometer_and_fiber_array(
    cell=right4x4_cell,
    inports=device_inports_2,
    gc_positions=grating_coupler_positions_2,
    electrode_length=1250,
    coupler_sep=0.5,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=9,
    wf_layer=100,
    electrode_layer=10,
    right_electrode_wf_layers=right_electrode_wf_layers,
    left_electrode_wf_layers=left_electrode_wf_layers,
    electrode_contact_pad_coordinates=[],
    second_x_middle=second_x_middle_2,
    delay_lines=False,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    last_interferometer=True,
    mm_wg_width=wg_Expwidth,
    mm_taper_length=30,
    contact_pad_dims=[250, 250],
    contact_pad_sep=160,
    connector_start_y=7180,
    connector_probe_pitch=175,
    wg_sep=wg_sep,
    electrode_wg_sep=electrode_wg_sep,
    gc_params={
            'width': 0.5,
            'full_opening_angle': np.deg2rad(180),
            'grating_period': 0.5,
            'grating_ff': 0.2,
            'n_gratings': 10,
            'ap_max_ff': 0.8,
            'n_ap_gratings': 65,
            'taper_length': 9.64
    }
)
right4x4_cell.add_to_layer(comment_layer, Text((4060, 540), 10,
                    'wg_width={:.2f}um\nexp_wg_width={:.2f}um\ngrat_period={:.2f}um'.format(wg_width,wg_Expwidth, 0.5),
                    alignment='center-top'))


# Adding global markers
global_markers_cell = Cell("Global markers")
build_global_markers(cell=global_markers_cell,
                     bounds=[-1015, -160, 10240, 8700],
                     marker_diameter=20,
                     marker_layer=2,
                     wg_layer=9,
                     marker_protection_layer=15)

# #######################################################
# ### ADDING TIME-BIN INTERFEROMETERS

timebin_cell = Cell('TimeBin_Interferometers')
#
x_adj = 6368.508 - 6330.45
y_adj = -130.242 + 150.239

OriginTests = (first_x_middle_2 + 3056 +x_adj, 80 + y_adj)

scan_list = np.linspace(4, 32, 4)
# dev_dist_list = (0, 310+(7212.45 - 7255.508), 410 + (8194.45-8199.45), 510 + (9276.45- 9281.45))
dev_dist_list = (0, 310+(7212.45 - 7178.508), 410 + (8194.45-8122.45), 510 + (9276.45- 9204.45))
dx_adj_list = (490.616-571.668 + (1.124- 0.374), 528.674 - 571.668 + (1.124- 0.374), 528.674 - 571.668 + (1.124- 0.374), 528.674 - 571.668 + (1.124- 0.374))
add_ylength_list = [6350,7750,7750,7750]


temp_pos = OriginTests
for x, num in enumerate(scan_list):
    temp_cell, x_max = TimeBin_BS('BS%i' % x, 2 * int(num / 2.), add_xlength=0., add_ylength=add_ylength_list[x], sep=5., r_curve=70.,
                                  coupler_sep=0.5, coupler_length=30, MZ_length=1250, electrodes_sep=1.1,
                                  return_xmax=True, dx_adj = dx_adj_list[x], exp_wg_width=wg_Expwidth, grating_coupler_period = std_coupler_params['grating_period'])
    timebin_cell.add_cell(temp_cell, origin=temp_pos + np.array((dev_dist_list[x], 0)), angle=np.pi)
    temp_pos = np.array((temp_pos[0] + x_max + dev_dist_list[x], OriginTests[1]))


# #######################################################
# ### ADDING DEMUX SETUP

demux_cell = Cell('Demux')

OriginTests = (1460, 1120)

coupler_sep = 0.5
coupler_length = 30 # old value
electrode_length = 1250

electrodes_sep = mod_params['electrode_sep']
row = 0

# MZIlen_list = [700, 975, 1250, 1515]
MZIlen_list = [1250]
x_dist = 930
y_dist = 1760
# y_steps = 275
y_steps = -265

for col, MZ_length in enumerate(MZIlen_list):
    temp_cell = Demux_active(coupler_sep, coupler_length, MZ_length, electrodes_sep,
                           id_to_alphanumeric(row, col), exp_wg_width=wg_Expwidth, grating_coupler_period = std_coupler_params['grating_period'])
    demux_cell.add_cell(temp_cell, origin=(OriginTests[0], OriginTests[1]), angle=np.pi)

# ### ADDING DEMUX WRITEFIELDS
wf_maxlength = 1040
dmx_wf_layer = 210

first_line_bounds = [1020, 960, 1650, 1270]
second_line_bounds = [1180, 1270, 1420, 3175]
turn_bounds = [650, 3175, 1420, 3290]
third_line_bounds = [650, 1270, 880, 3175]
fourth_line_bounds = [370, 960, 1010, 1270]

dmx_region_marker = bounds_to_polygon(first_line_bounds)

_, dmx_region_marker = wf_line_from_bounds(
    bounds=first_line_bounds,
    cell=demux_cell,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)
_, dmx_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=second_line_bounds,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)
_, dmx_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=turn_bounds,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)
_, dmx_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=third_line_bounds,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)
_, dmx_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=fourth_line_bounds,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)

connect_line = [
    fourth_line_bounds[0],
    left_in_bounds[3],
    fourth_line_bounds[0]+20,
    fourth_line_bounds[1]
]

_, dmx_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=connect_line,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)

## Very brute-force
left_in_bounds = [left_in_bounds[0]+314,
                  left_in_bounds[1],
                  left_in_bounds[2],
                  left_in_bounds[3]]

_, dmx_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=left_in_bounds,
    region_marker=dmx_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=dmx_wf_layer
)

demux_cell.add_to_layer(dmx_wf_layer+1, dmx_region_marker)

### ELECTRODES

el1_first_line_bounds = [330,1500, 870, 3241]
el1_second_line_bounds = [330, 3241, 870, 4010]
el1_third_line_bounds = [870,3241,1660,4010]
el1_fourth_line_bounds = [1200,1520,1590,3241]


el1_wf_layer = 220

el1_region_marker = bounds_to_polygon(el1_first_line_bounds)

_, el1_region_marker = wf_line_from_bounds(
    bounds=el1_first_line_bounds,
    cell=demux_cell,
    region_marker=el1_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=el1_wf_layer
)
_, el1_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=el1_second_line_bounds,
    region_marker=el1_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=el1_wf_layer
)
_, el1_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=el1_third_line_bounds,
    region_marker=el1_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=el1_wf_layer
)
_, el1_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=el1_fourth_line_bounds,
    region_marker=el1_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=el1_wf_layer
)

connect_line = [
    el1_first_line_bounds[0],
    left_in_bounds[3],
    el1_first_line_bounds[0]+20,
    el1_first_line_bounds[1]
]
_, el1_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=connect_line,
    region_marker=el1_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=el1_wf_layer
)

_, el1_region_marker = wf_line_from_bounds(
    cell=demux_cell,
    bounds=left_in_bounds,
    region_marker=el1_region_marker,
    wf_maxlength=wf_maxlength,
    wf_layer=el1_wf_layer
)

demux_cell.add_to_layer(el1_wf_layer+1, el1_region_marker)


##############################################################
### Add all test structures

tests_cell = Cell('Test_Structures')

## ADDING DELAY_LINES TEST STRUCTURES
init_spirals_pos = np.array((390, 4190))

scan_list = np.linspace(4, 40, 4)
dev_dist_list = (0, 310, 310, 360)
grating_angle = np.deg2rad(180)

temp_pos = init_spirals_pos
for x, num in enumerate(scan_list):
    temp_cell, x_max = RectangularSpiral(id_to_alphanumeric(x, 0), 2 * int(num / 2.), add_xlength=0, add_ylength=1100,
                                         sep=5, r_curve=80, return_xmax=True, grating_angle=grating_angle, exp_wg_width=wg_Expwidth, grating_coupler_period = std_coupler_params['grating_period'])
    tests_cell.add_cell(temp_cell, origin=temp_pos + np.array((dev_dist_list[x], 0)))
    temp_pos = np.array((temp_pos[0] + x_max + dev_dist_list[x], init_spirals_pos[1]))

# ### ADDING DIRECTIONAL COUPLERS TEST STRUCTURES - STANDARD ONES

# print(dL)
spacingy = np.array((0, 500))

# right
OriginDCTests = np.array((3200, 440))
col = 0
gap = 0.5
NumTests_y = 3
for row, length in enumerate(np.linspace(28.5, 31.5, NumTests_y)):
    temp_cell = DirectionalCouplersTest_standard(str(col) + ',' + str(row), gap, length)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy, angle=np.pi)


# centre
OriginDCTests = np.array((2490, 940))
col = 1
gap = 0.48
NumTests_y = 2
for row, length in enumerate(np.linspace(28.5, 30, NumTests_y)):
    temp_cell = DirectionalCouplersTest_standard(str(col) + ',' + str(row), gap, length)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy, angle=np.pi)

# left
OriginDCTests = np.array((-325, 952))
col = 2
gap = 0.52
NumTests_y = 2
for row, length in enumerate(np.linspace(28.5, 30, NumTests_y)):
    temp_cell = DirectionalCouplersTest_standard(str(col) + ',' + str(row), gap, length)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy, angle=np.pi)

# ### ADDING GRATING COUPLERS TEST STRUCTURES

# left side - bottom
OriginTests = np.array((-280, 1850))
spacingy = np.array((0, 170))
col = 0
ff = 0.15
NumTests_y = 7
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# left side - top
OriginTests = np.array((-280, 3180))
spacingy = np.array((0, 170))
col = 1
ff = 0.2
NumTests_y = 8
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# central side - right
OriginTests = np.array((1772, 120))
spacingy = np.array((0, 170))
col = 2
ff = 0.25
NumTests_y = 5
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff, grating_angle = np.deg2rad(120))
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side
OriginTests = np.array((2220, 3380))
spacingy = np.array((0, 170))
col = 3
ff = 0.3
NumTests_y = 7
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)


# ### ADDING RING RESONATORS - RIGHT 4X4 SPACES

# left side - centre
spacingy = np.array((0, 220))
OriginTests = np.array((first_x_middle_2 - 300, 1940))  # (3200, 1850) for probe pitch 150
col = 0
ring_r = 30
NumTests_y = 5
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# left side - top
spacingy = np.array((0, 240))
OriginTests = np.array((first_x_middle_2 - 300, 3300))  # (3200, 3200) for probe pitch 150
col = 1
ring_r = 40
NumTests_y = 5
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# # right side - centre left
spacingy = np.array((0, 322))
OriginTests = np.array((first_x_middle_2 + 1285 -2500, 2090))  # (5700, 1950) for probe pitch 150
col = 2
ring_r = 50
NumTests_y = 4
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# # right side - centre left
spacingy = np.array((0, 322))
OriginTests = np.array((first_x_middle_2 + 2185, 2040))  # (5700, 1950) for probe pitch 150
col = 3
ring_r = 60
NumTests_y = 4
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - centre right
spacingy = np.array((0, 322))
OriginTests = np.array((first_x_middle_2 + 2430, 2040))  # (5700, 1950) for probe pitch 150
col = 4
ring_r = 55
NumTests_y = 4
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - bottom left
spacingy = np.array((0, 302))
OriginTests = np.array((first_x_middle_2 + 2220, 580))  # (5920, 1960) for probe pitch 150
col = 5
ring_r = 70
NumTests_y = 4
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - bottom right
spacingy = np.array((0, 302))
OriginTests = np.array((first_x_middle_2 + 2500, 580))  # (5920, 1960) for probe pitch 150
col = 6
ring_r = 65
NumTests_y = 4
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - top
spacingy = np.array((0, 362))
OriginTests = np.array((first_x_middle_2 + 2200, 3570))  # (5700, 3450) for probe pitch 150
col = 7
ring_r = 80
NumTests_y = 3
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r, wg_Expwidth)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

### ADDING MZI TESTS
OriginTests = (first_x_middle_2 + 800, 1650-120)

coupler_sep = 0.5
coupler_length = 30 # old value
# coupler_length = 14.5  # New value from francesco
electrode_length = 1250

electrodes_sep = mod_params['electrode_sep']
row = 0

MZIlen_list = [650, 925, 1250, 1465]
x_dist = 930
y_dist = 1760
y_steps = 275

for col, MZ_length in enumerate(MZIlen_list):
    temp_cell = MZI_active(coupler_sep, coupler_length, MZ_length, electrodes_sep,
                           id_to_alphanumeric(row, col), exp_wg_width=wg_Expwidth, grating_coupler_period = std_coupler_params['grating_period'])
    tests_cell.add_cell(temp_cell,
                        origin=(OriginTests[0] + (col % 2) * x_dist , OriginTests[1] + int(col / 2.) * y_dist +
                                (col % 2) * int(col / 2.) * y_steps), angle=np.pi)




##############################################################
### ADD ALL CELLS AND EXPORT FULL GDS
global_cell = Cell('LiNb01_MunsterNBI')

global_cell.add_cell(left4x4_cell, origin=(0, 0))
global_cell.add_cell(right4x4_cell, origin=(0, 0))
global_cell.add_cell(timebin_cell, origin=(0, 0))
global_cell.add_cell(demux_cell, origin=(0, 0))
global_cell.add_cell(tests_cell, origin=(0, 0))
global_cell.add_cell(global_markers_cell, origin=(0, 0))


global_cell.save('LiNb01_MunsterNBI_v5_500nm.gds')
