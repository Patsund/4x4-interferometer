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
from parts import interferometer_and_fiber_array
from TestStructures.tests_DC_Grating_Delays import RectangularSpiral, DirectionalCouplersTest, Efficiency_Grating



cell = Cell('4x4 interferometer')

#Layers specs
wg_layer = 9
box_layer = 100
comment_layer = 11
marker_layer_1 = 3
marker_layer_2 = 4
marker_protection_layer = 15

wg_width = 0.5
electrode_wg_sep = 62.5
second_x_middle = 8 * 2 * 127

inport_0 = Port((-electrode_wg_sep/2 - 25, 0), np.pi/2, wg_width)
inport_1 = Port((-electrode_wg_sep/2, 0), np.pi/2, wg_width)
inport_2 = Port((electrode_wg_sep/2, 0), np.pi/2, wg_width)
inport_3 = Port((electrode_wg_sep/2 + 25, 0), np.pi/2, wg_width)
device_inports = [inport_0, inport_1, inport_2, inport_3]

start_x = inport_0.origin[0] + 3500
device_inports_2 = [None, None, None, None]
device_inports_2[:2] = [
    Port((start_x + idx*25, 0), np.pi/2, wg_width)
    for idx in range(2)
]
start_x += 25 + electrode_wg_sep
first_x_middle_2 = start_x - electrode_wg_sep/2
second_x_middle_2 = first_x_middle_2 + second_x_middle
device_inports_2[2:] = [
    Port((start_x + idx*25, 0), np.pi/2, wg_width)
    for idx in range(2)
]

# fiber_array
gc_pitch = 127
wg_sep = 25
min_radius = 50
gc_starting_point = second_x_middle/2 - gc_pitch * 3.5
gc_starting_point_2 = (
    0.5*(first_x_middle_2 + second_x_middle_2)
    - gc_pitch * 3.5
)
gc_y = 100
grating_coupler_positions = [(gc_starting_point+idx*gc_pitch, gc_y)
                             for idx in range(8)]
grating_coupler_positions_2 = [(gc_starting_point_2+idx*gc_pitch, gc_y)
                               for idx in range(8)]

# Adding Left 4x4 interferometer
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

# Adding Right 4x4 interferometer
_ = interferometer_and_fiber_array(
    cell=cell,
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

### ADDING TEST STRUCTURES
devices = []

### Adding Delay-Lines test structures
init_spirals_pos = np.array((405, 980))

scan_list = np.linspace(4, 36, 5)
dev_dist_list = (0, 234, 234, 262, 302)

temp_pos = init_spirals_pos
for x, num in enumerate(scan_list):
    this_dev = RectangularSpiral(temp_pos + np.array((dev_dist_list[x], 0)), id_to_alphanumeric(x, 0),
                                 2 * int(num / 2.),
                                 add_xlength=0, add_ylength=4350, sep=5, r_curve=45)
    devices += [this_dev]
    this_box = (this_dev.get_box())[0].get_shapely_object()
    temp_pos = np.array((max(this_box.exterior.coords.xy[0]), init_spirals_pos[1]))


### Adding directional couplers test structures

def FSR_to_L(FSR_GHz, n_g=2.1):
    # CHECK THAT THE REFRACTIVE INDEX IS CORRECT!!!!
    c = 2.998 * 10 ** 8
    return c / n_g / FSR_GHz * 10 ** 6 * 10 ** -9

dL = FSR_to_L(160)

print(dL)
spacingy = np.array((0, 212))

#left side
OriginDCTests = np.array((-500, -5))
col=0
gap=0.3
NumTests_y = 8
for row, length in enumerate(np.linspace(15, 35, NumTests_y)):
    devices.append(
        DirectionalCouplersTest(OriginDCTests + row * spacingy, str(col) + ',' + str(row), gap,
                                length, dL, r_curve=35))
#right-bottom side
OriginDCTests = np.array((2220, 445))
col=1
gap=0.4
NumTests_y = 6
for row, length in enumerate(np.linspace(20, 35, NumTests_y)):
    devices.append(
        DirectionalCouplersTest(OriginDCTests + row * spacingy, str(col) + ',' + str(row), gap,
                                length, dL, r_curve=35))

#right-top side
OriginDCTests = np.array((2220, 1910))
col=2
gap=0.5
NumTests_y = 6
for row, length in enumerate(np.linspace(20, 35, NumTests_y)):
    devices.append(
        DirectionalCouplersTest(OriginDCTests + row * spacingy, str(col) + ',' + str(row), gap,
                                length, dL, r_curve=35))


### Adding grating couplers test structures
spacingy = np.array((0, 215))

#left side - bottom
OriginTests = np.array((-280, 1870))
col = 0
ff = 0.25
NumTests_y = 5
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    devices.append(
        Efficiency_Grating(OriginTests + row * spacingy, id_to_alphanumeric(row, col+1), period, ff)
        )

#left side - top
OriginTests = np.array((-280, 3120))
col = 1
ff = 0.3
NumTests_y = 6
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    devices.append(
        Efficiency_Grating(OriginTests + row * spacingy, id_to_alphanumeric(row, col+1), period, ff)
        )

#central side - left
OriginTests = np.array((1555, 150))
col = 2
ff = 0.2
NumTests_y = 4
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    devices.append(
        Efficiency_Grating(OriginTests + row * spacingy, id_to_alphanumeric(row, col+1), period, ff)
        )

#central side - right
OriginTests = np.array((1757, 150))
col = 3
ff = 0.4
NumTests_y = 4
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    devices.append(
        Efficiency_Grating(OriginTests + row * spacingy, id_to_alphanumeric(row, col+1), period, ff)
        )

#right side
OriginTests = np.array((2220, 3350))
col = 4
ff = 0.35
NumTests_y = 5
for row, period in enumerate(np.linspace(0.4, 0.52, NumTests_y)):
    devices.append(
        Efficiency_Grating(OriginTests + row * spacingy, id_to_alphanumeric(row, col+1), period, ff)
        )



### Add all test structures
for i, device in enumerate(devices):
    cell.add_to_layer(wg_layer, device.get_shapely_object())
    device_boxes = device.get_box()
    if isinstance(device_boxes, list):
        for this_box in device.get_box():
            cell.add_to_layer(box_layer, this_box)
    else:
        cell.add_to_layer(box_layer, device_boxes)
    cell.add_to_layer(comment_layer, device.get_comments().get_shapely_object())
    cell.add_to_layer(marker_layer_1, device.get_markers())
    cell.add_to_layer(marker_layer_2, device.get_markers_copy())
    cell.add_to_layer(marker_protection_layer, device.get_markers_protection())



cell.save('4x4_two_devices.gds')
