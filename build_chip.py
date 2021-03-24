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

from tech.LiNb01 import *

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
min_radius = 50
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

_ = interferometer_and_fiber_array(
    cell=left4x4_cell,
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
    connector_start_y=7180,
    connector_probe_pitch=175,
    wg_sep=wg_sep,
    electrode_wg_sep=electrode_wg_sep
)

# Adding Right 4x4 interferometer
right4x4_cell = Cell('4x4 interferometer_Right')
_ = interferometer_and_fiber_array(
    cell=right4x4_cell,
    inports=device_inports_2,
    gc_positions=grating_coupler_positions_2,
    electrode_length=1250,
    coupler_sep=0.4,
    coupler_length=30,
    sm_wg_width=wg_width,
    wg_layer=9,
    wf_layer=100,
    electrode_layer=10,
    electrode_wf_layer=102,
    electrode_contact_pad_coordinates=[],
    second_x_middle=second_x_middle_2,
    delay_lines=False,
    total_bounds=(np.inf, np.inf, -np.inf, -np.inf),
    mm_wg_width=1,
    mm_taper_length=30,
    contact_pad_dims=[250, 250],
    contact_pad_sep=160,
    connector_start_y=7180,
    connector_probe_pitch=175,
    wg_sep=wg_sep,
    electrode_wg_sep=electrode_wg_sep
)

# Adding global markers
global_markers_cell = Cell("Global markers")
build_global_markers(cell=global_markers_cell,
                     bounds=[-815, -160, 10240, 8800],
                     marker_diameter=20,
                     marker_layer=2,
                     wg_layer=9,
                     marker_protection_layer=15)

#######################################################
### ADDING TIME-BIN INTERFEROMETERS

timebin_cell = Cell('TimeBin_Interferometers')

OriginTests = (first_x_middle_2 + 3056, 80)

scan_list = np.linspace(4, 32, 4)
dev_dist_list = (0, 360, 460, 560)

temp_pos = OriginTests
for x, num in enumerate(scan_list):
    temp_cell, x_max = TimeBin_BS('BS%i' % x, 2 * int(num / 2.), add_xlength=0., add_ylength=8000., sep=5., r_curve=50.,
                                  coupler_sep=0.4, coupler_length=30, MZ_length=1250, electrodes_sep=1.1,
                                  return_xmax=True)
    timebin_cell.add_cell(temp_cell, origin=temp_pos + np.array((dev_dist_list[x], 0)), angle=np.pi)
    temp_pos = np.array((temp_pos[0] + x_max + dev_dist_list[x], OriginTests[1]))

##############################################################
### Add all test structures

tests_cell = Cell('Test_Structures')

### ADDING DELAY_LINES TEST STRUCTURES
init_spirals_pos = np.array((405, 980))

scan_list = np.linspace(4, 36, 5)
dev_dist_list = (0, 234, 234, 262, 302)

temp_pos = init_spirals_pos
for x, num in enumerate(scan_list):
    temp_cell, x_max = RectangularSpiral(id_to_alphanumeric(x, 0), 2 * int(num / 2.), add_xlength=0, add_ylength=4350,
                                         sep=5, r_curve=45, return_xmax=True)
    tests_cell.add_cell(temp_cell, origin=temp_pos + np.array((dev_dist_list[x], 0)))
    temp_pos = np.array((temp_pos[0] + x_max + dev_dist_list[x], init_spirals_pos[1]))


### ADDING DIRECTIONAL COUPLERS TEST STRUCTURES

def FSR_to_L(FSR_GHz, n_g=2.1):
    # CHECK THAT THE REFRACTIVE INDEX IS CORRECT!!!!
    c = 2.998 * 10 ** 8
    return c / n_g / FSR_GHz * 10 ** 6 * 10 ** -9


dL = FSR_to_L(160)

# print(dL)
spacingy = np.array((0, 212))

# left side
OriginDCTests = np.array((-500, -5))
col = 0
gap = 0.3
NumTests_y = 8
for row, length in enumerate(np.linspace(15, 35, NumTests_y)):
    temp_cell = DirectionalCouplersTest(str(col) + ',' + str(row), gap, length, dL, r_curve=35)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy)

# right-bottom side
OriginDCTests = np.array((2220, 445))
col = 1
gap = 0.4
NumTests_y = 6
for row, length in enumerate(np.linspace(20, 35, NumTests_y)):
    temp_cell = DirectionalCouplersTest(str(col) + ',' + str(row), gap, length, dL, r_curve=35)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy)

# right-top side
OriginDCTests = np.array((2220, 1910))
col = 2
gap = 0.5
NumTests_y = 6
for row, length in enumerate(np.linspace(20, 35, NumTests_y)):
    temp_cell = DirectionalCouplersTest(str(col) + ',' + str(row), gap, length, dL, r_curve=35)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy)



### ADDING DIRECTIONAL COUPLERS TEST STRUCTURES - STANDARD ONES

# print(dL)
spacingy = np.array((0, 500))

# centre
OriginDCTests = np.array((3200, 400))
col = 0
gap = 0.4
NumTests_y = 3
for row, length in enumerate(np.linspace(25, 35, NumTests_y)):
    temp_cell = DirectionalCouplersTest_standard(str(col) + ',' + str(row), gap, length)
    tests_cell.add_cell(temp_cell, origin=OriginDCTests + row * spacingy, angle=np.pi)


### ADDING GRATING COUPLERS TEST STRUCTURES
spacingy = np.array((0, 100))

# left side - bottom
OriginTests = np.array((-280, 1760))
col = 0
ff = 0.25
NumTests_y = 11
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# left side - top
OriginTests = np.array((-280, 3070))
col = 1
ff = 0.3
NumTests_y = 12
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# central side - right
OriginTests = np.array((1757, 60))
col = 2
ff = 0.4
NumTests_y = 9
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side
OriginTests = np.array((2220, 3280))
col = 3
ff = 0.35
NumTests_y = 12
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)


### ADDING GRATING COUPLERS TEST STRUCTURES - FINE SCAN BETWEEN LONG STRUCTURES
spacingy = np.array((0, 100))

# left - bottom
OriginTests = np.array((first_x_middle_2 + 3485, 3000))
col = 10
ff = 0.25
NumTests_y = 27
for row, period in enumerate(np.linspace(0.32, 0.6, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# left - top
OriginTests = np.array((first_x_middle_2 + 3485, 6000))
col = 11
ff = 0.275
NumTests_y = 27
for row, period in enumerate(np.linspace(0.32, 0.6, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# centre - bottom
OriginTests = np.array((first_x_middle_2 + 4420, 3000))
col = 12
ff = 0.30
NumTests_y = 27
for row, period in enumerate(np.linspace(0.32, 0.6, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# centre - top
OriginTests = np.array((first_x_middle_2 + 4420, 6000))
col = 13
ff = 0.325
NumTests_y = 27
for row, period in enumerate(np.linspace(0.32, 0.6, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right - bottom
OriginTests = np.array((first_x_middle_2 + 5450, 3000))
col = 14
ff = 0.35
NumTests_y = 27
for row, period in enumerate(np.linspace(0.32, 0.6, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right - top
OriginTests = np.array((first_x_middle_2 + 5450, 6000))
col = 15
ff = 0.37
NumTests_y = 27
for row, period in enumerate(np.linspace(0.32, 0.6, NumTests_y)):
    temp_cell = Efficiency_Grating(id_to_alphanumeric(row, col + 1), period, ff)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)


### ADDING RING RESONATORS - RIGHT 4X4 SPACES

# left side - centre
spacingy = np.array((0, 160))
OriginTests = np.array((first_x_middle_2 - 300, 1850))  # (3200, 1850) for probe pitch 150
col = 0
ring_r = 30
NumTests_y = 7
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# left side - top
spacingy = np.array((0, 182))
OriginTests = np.array((first_x_middle_2 - 300, 3200))  # (3200, 3200) for probe pitch 150
col = 1
ring_r = 40
NumTests_y = 7
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# central side
spacingy = np.array((0, 202))
OriginTests = np.array((first_x_middle_2 + 1770, 180))  # (5220, 180) for probe pitch 150
col = 2
ring_r = 50
NumTests_y = 7
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - centre left
spacingy = np.array((0, 222))
OriginTests = np.array((first_x_middle_2 + 2200, 1950))  # (5700, 1950) for probe pitch 150
col = 3
ring_r = 60
NumTests_y = 6
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - centre right
spacingy = np.array((0, 222))
OriginTests = np.array((first_x_middle_2 + 2430, 1950))  # (5700, 1950) for probe pitch 150
col = 4
ring_r = 55
NumTests_y = 6
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - bottom left
spacingy = np.array((0, 242))
OriginTests = np.array((first_x_middle_2 + 2220, 540))  # (5920, 1960) for probe pitch 150
col = 5
ring_r = 70
NumTests_y = 5
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - bottom right
spacingy = np.array((0, 242))
OriginTests = np.array((first_x_middle_2 + 2500, 540))  # (5920, 1960) for probe pitch 150
col = 6
ring_r = 65
NumTests_y = 5
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

# right side - top
spacingy = np.array((0, 262))
OriginTests = np.array((first_x_middle_2 + 2200, 3500))  # (5700, 3450) for probe pitch 150
col = 7
ring_r = 80
NumTests_y = 4
for row, gap in enumerate(np.linspace(0.2, 0.8, NumTests_y)):
    temp_cell = Ring_Test(id_to_alphanumeric(row, col), gap, ring_r)
    tests_cell.add_cell(temp_cell, origin=OriginTests + row * spacingy)

### ADDING MZI TESTS
OriginTests = (first_x_middle_2 + 775, 1610)

coupler_sep = 0.4
coupler_length = 30
electrode_length = 1250

electrodes_sep = mod_params['electrode_sep']
row = 0

MZIlen_list = [700, 975, 1250, 1515]
x_dist = 930
y_dist = 1760
y_steps = 275
for col, MZ_length in enumerate(MZIlen_list):
    temp_cell = MZI_active(coupler_sep, coupler_length, MZ_length, electrodes_sep,
                           id_to_alphanumeric(row, col))
    tests_cell.add_cell(temp_cell,
                        origin=(OriginTests[0] + (col % 2) * x_dist, OriginTests[1] + int(col / 2.) * y_dist +
                                (col % 2) * int(col / 2.) * y_steps), angle=np.pi)

# ##############################################################
# ### ADD ALL CELLS AND EXPORT FULL GDS
global_cell = Cell('LiNb01_MunsterNBI')

global_cell.add_cell(left4x4_cell, origin=(0, 0))
global_cell.add_cell(right4x4_cell, origin=(0, 0))
global_cell.add_cell(timebin_cell, origin=(0, 0))
global_cell.add_cell(tests_cell, origin=(0, 0))
global_cell.add_cell(global_markers_cell, origin=(0, 0))

global_cell.save('LiNb01_MunsterNBI_v0.gds')
