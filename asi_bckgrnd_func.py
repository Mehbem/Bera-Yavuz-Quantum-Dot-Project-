import zwoasi as asi
import os
import sys 
import cv2
import numpy as np
import matplotlib.pyplot as plt
import qd_data_folder_creation

def snap_background(asi_initialization_return):
    # Fetch today's folder string 
    date_string = str(qd_data_folder_creation.date_string())
    date_string_test = date_string+"_Test"

    # Enter directory where the images are being saved
    qd_data_directory = r"c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data"
    ASI_raw_directory = str(date_string_test)+"\Spectrometer_ASI"
    qd_data_ASI_raw_directory = os.path.join(qd_data_directory, ASI_raw_directory)
    os.chdir(qd_data_ASI_raw_directory) 

    filename = 'background_1800mm_grating_'+date_string

    camera = asi_initialization_return

    for i in range(3):
        camera_img_filename = filename + '_' + str(i + 1) + '.png'
        camera.set_image_type(asi.ASI_IMG_RAW8)
        camera.capture(filename=camera_img_filename)
        print('Saved to %s' % camera_img_filename)

