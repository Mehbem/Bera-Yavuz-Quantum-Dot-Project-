import zwoasi as asi
import os
import sys 
import cv2
import numpy as np
import matplotlib.pyplot as plt
import mplfig
import qd_data_folder_creation

def init_camera():

    # File downloaded from the ASI SDK website
    asi.init('ASICamera2.lib')

    # Camera settings
    default_exposure = 2000000
    exposure_setting = default_exposure #int(0.0078E6)

    # Initialize camera
    num_cameras = asi.get_num_cameras()
    if num_cameras == 0:
        raise ValueError('No cameras found')

    camera_id = 0  # use first camera from list
    cameras_found = asi.list_cameras()
    print(cameras_found)
    camera = asi.Camera(camera_id)
    camera_info = camera.get_camera_property()
    # print(camera_info)

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

    # Begin image taking 
    print('Enabling stills mode')
    try:
        # Force any single exposure to be halted
        camera.stop_video_capture()
        camera.stop_exposure()
    except (KeyboardInterrupt, SystemExit):
        raise

    asi_initialization_return = camera

    return asi_initialization_return

def snap_image(asi_initialization_return, filename, dir_bckgrnd, dir_save_ASI, dir_save_plots):

    # Fetch today's folder string 
    date_string = str(qd_data_folder_creation.date_string())
    date_string_test = date_string+"_Test"

    # Parameters for different gratings
    row_index_1800mm_grating = 2000
    row_index_1200mm_grating = 2495
    row_index_150mm_grating = 1795
    row_index_setting = row_index_1200mm_grating

    # Read background images and obtain the background count arrays
    os.chdir(dir_bckgrnd) 
    bckgrnd_avg_array_filename = 'background_1200mm_grating_avg_array_'+ date_string+'.npy'
    bckgrnd_avg = np.load(bckgrnd_avg_array_filename)

    # Save new background 
    os.chdir(dir_save_ASI) 
    camera = asi_initialization_return

    # print('Capturing a single 8-bit mono image')
    camera_img_filename = filename + '.png'
    camera.set_image_type(asi.ASI_IMG_RAW8)
    camera.capture(filename=camera_img_filename)
    #print('Saved to %s' % camera_img_filename)

    # OpenCV checks, read and save image
    img = cv2.imread(camera_img_filename, cv2.IMREAD_GRAYSCALE)
    y_index_range = np.arange(-400, 400, 1)
    spectrum_sum = np.zeros(len(img[1, :]))

    #plt.imshow(img)
    #plt.show()

    wvlngth_start = 889.1714 
    wvlngth_end =  900.8576	

    dummy_wvlength = np.linspace(wvlngth_start, wvlngth_end, len(img[1, :]))

    for j in range(len(y_index_range)):
        spectrum_sum = spectrum_sum + img[row_index_setting + y_index_range[j], :]

    spectrum_sum = spectrum_sum - bckgrnd_avg 
    spectrum_sum = spectrum_sum/2
    
    # Save plot .png sile of the spectrum plot 
    os.chdir(dir_save_plots) 

    fig, ax = plt.subplots()
    ax.set_box_aspect(0.5)
    ax.plot(dummy_wvlength, spectrum_sum)
    ax.set_title('QD Spectrum '+filename)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Arb. Counts')
    plt.gcf().set_size_inches(12, 8)
    #plt.show()
    plt.savefig(filename+'_wvl_graph'+'.png')
    plt.clf()

    cv2.destroyAllWindows()

def snap_background(asi_initialization_return, dir_bckgrnd):
    # Fetch today's folder string 
    date_string = str(qd_data_folder_creation.date_string())
    date_string_test = date_string+"_Test"

    # Change directory for saving backgrounds 
    os.chdir(dir_bckgrnd) 
    
    bckgrnd_filename = 'background_1200mm_grating_'+date_string
    camera = asi_initialization_return

    for i in range(3):
        camera_img_filename = bckgrnd_filename + '_' + str(i + 1) + '.png'
        camera.set_image_type(asi.ASI_IMG_RAW8)
        camera.capture(filename=camera_img_filename)
        print('Saved to %s' % camera_img_filename)

    # Extract background data as an array 
    bckgrnd_img_filename_1 = bckgrnd_filename+'_1.png'
    bckgrnd_img_filename_2 = bckgrnd_filename+'_2.png'
    bckgrnd_img_filename_3 = bckgrnd_filename+'_3.png'

    bckgrnd_img_1 = cv2.imread(bckgrnd_img_filename_1, cv2.IMREAD_GRAYSCALE)
    bckgrnd_img_2 = cv2.imread(bckgrnd_img_filename_2, cv2.IMREAD_GRAYSCALE)
    bckgrnd_img_3 = cv2.imread(bckgrnd_img_filename_3, cv2.IMREAD_GRAYSCALE)

    bckgrnd_sum_1 = np.zeros(len(bckgrnd_img_1[1, :]))
    bckgrnd_sum_2 = np.zeros(len(bckgrnd_img_1[1, :]))
    bckgrnd_sum_3 = np.zeros(len(bckgrnd_img_1[1, :]))

    # Parameters for different gratings
    row_index_1800mm_grating = 2000
    row_index_1200mm_grating = 2495
    row_index_150mm_grating = 1795
    row_index_setting = row_index_1200mm_grating
    
    y_index_range = np.arange(-400, 400, 1)

    for j in range(len(y_index_range)):
        bckgrnd_sum_1 = bckgrnd_sum_1 + bckgrnd_img_1[row_index_setting + y_index_range[j], :]
        bckgrnd_sum_2 = bckgrnd_sum_2 + bckgrnd_img_2[row_index_setting + y_index_range[j], :]
        bckgrnd_sum_3 = bckgrnd_sum_3 + bckgrnd_img_3[row_index_setting + y_index_range[j], :]

    bckgrnd_avg = (bckgrnd_sum_1 + bckgrnd_sum_2 + bckgrnd_sum_3)/3

    bckgrnd_avg_array_filename = 'background_1200mm_grating_avg_array_'+ date_string+'.npy'
    np.save(bckgrnd_avg_array_filename, bckgrnd_avg)


#asi_initialization_return = init_camera() 
#snap_background(asi_initialization_return)
#snap_image(asi_initialization_return, 'TEST_6s_QD_Spec_Plot_[78 6]')