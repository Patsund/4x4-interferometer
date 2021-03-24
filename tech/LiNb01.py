import numpy as np

### Waveguides
wg_width = 0.5
wg_Expwidth = 1.0
l_Exptaper = 40

### Bend radii
bend_r = 50
euler_to_bend_coeff = 1.8703879865

### Markers
marker_dims = 20

### Box Sizes
max_box_size = (1040, 1040)

### Layers
wg_layer = 9
electrode_layer = 10
wg_wf_layer = 100
wg_reg_layer = 101
electrode_wf_layer = 102
electrode_reg_layer = 103
comment_layer = 6
marker_layer_1 = 3
marker_layer_2 = 4
marker_protection_layer = 15

### Grating Couplers
std_coupler_params = {
    'width': 0.5,
    'full_opening_angle': np.deg2rad(90),
    'grating_period': 0.46,
    'grating_ff': 0.3,
    'n_gratings': 10,
    'ap_max_ff': 0.8,
    'n_ap_gratings': 55,
    'taper_length': 12
}

### Modulators
mod_params = {
    'wg_sep': 25,
    'mm_taper_length': 0,
    'sine_s_x': 60,
    'wf_maxlength': 1040,
    'wf_leeways': (10, 10),
    # For electrode function
    'electrode_width': 25,
    'electrode_sep': 1.1,
    'crossing_width': 10,
    'electrode_taper_leeway': 5,
    'electrode_taper_length': 30,
    'electrode_sep_y': 40,
    'connector_probe_dims': (80, 750),
    'connector_probe_pitch': 175,
    'contact_pad_dims': (250, 600),
    'contact_pad_sep': 175,
    'electrode_pad_seps': (30, 30),
    'pad_sep': 300,
}


### Others
opt_space = 127
box_dist = 5

