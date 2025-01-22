import zwoasi as asi
import os
import sys 
import cv2
import numpy as np
import matplotlib.pyplot as plt
import qd_data_folder_creation

# Fetch today's folder string 
date_string = str(qd_data_folder_creation.date_string())
date_string_test = date_string+"_Test"

# Enter directory where the images are being saved
qd_data_directory = r"c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data"
ASI_raw_directory = str(date_string_test)+"\Spectrometer_ASI"
qd_data_ASI_raw_directory = os.path.join(qd_data_directory, ASI_raw_directory)
os.chdir(qd_data_ASI_raw_directory) 

filename = 'background_1800mm_grating_'+date_string

# Camera settings
default_exposure = 2000000
exposure_setting = default_exposure #int(0.0078E6)

# File downloaded from the ASI SDK website
asi.init('ASICamera2.lib')

num_cameras = asi.get_num_cameras()
if num_cameras == 0:
    raise ValueError('No cameras found')

camera_id = 0  # use first camera from list
cameras_found = asi.list_cameras()
print(cameras_found)
camera = asi.Camera(camera_id)
camera_info = camera.get_camera_property()
#print(camera_info)

# Use minimum USB bandwidth permitted
camera.set_control_value(asi.ASI_BANDWIDTHOVERLOAD, camera.get_controls()['BandWidth']['MinValue'])

# Set some sensible defaults. They will need adjusting depending upon
# the sensitivity, lens and lighting conditions used.
camera.disable_dark_subtract()

camera.set_control_value(asi.ASI_GAIN, 300)
camera.set_control_value(asi.ASI_EXPOSURE, exposure_setting) # microseconds
camera.set_control_value(asi.ASI_WB_B, 99)
camera.set_control_value(asi.ASI_WB_R, 75)
camera.set_control_value(asi.ASI_GAMMA, 50)
camera.set_control_value(asi.ASI_BRIGHTNESS, 50)
camera.set_control_value(asi.ASI_FLIP, 0)

print('Enabling stills mode')
try:
    # Force any single exposure to be halted
    camera.stop_video_capture()
    camera.stop_exposure()
except (KeyboardInterrupt, SystemExit):
    raise

# print('Capturing a single 8-bit mono image')
for i in range(3):
    camera_img_filename = filename + '_' + str(i + 1) + '.png'
    camera.set_image_type(asi.ASI_IMG_RAW8)
    camera.capture(filename=camera_img_filename)
    print('Saved to %s' % camera_img_filename)