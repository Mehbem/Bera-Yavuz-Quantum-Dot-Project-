import numpy as np 

def read_spectrum(folder_name, camera_img_filename, bckgrnd_filename):

    # Background filename
    bckgrnd_img_filename_1 = bckgrnd_filename+'_1.png'
    bckgrnd_img_filename_2 = bckgrnd_filename+'_2.png'
    bckgrnd_img_filename_3 = bckgrnd_filename+'_3.png'

    # OpenCV checks, read and save image
    img = cv2.imread(camera_img_filename, cv2.IMREAD_GRAYSCALE)
    y_index_range = np.arange(-300, 300, 1)
    spectrum_sum = np.zeros(len(img[1, :]))

    # plt.imshow(img)
    # plt.show()

    bckgrnd_img_1 = cv2.imread(bckgrnd_img_filename_1, cv2.IMREAD_GRAYSCALE)
    bckgrnd_img_2 = cv2.imread(bckgrnd_img_filename_2, cv2.IMREAD_GRAYSCALE)
    bckgrnd_img_3 = cv2.imread(bckgrnd_img_filename_3, cv2.IMREAD_GRAYSCALE)

    bckgrnd_sum_1 = np.zeros(len(img[1, :]))
    bckgrnd_sum_2 = np.zeros(len(img[1, :]))
    bckgrnd_sum_3 = np.zeros(len(img[1, :]))

    wvlngth_start = 891.9909 #892.0069 #889.452
    wvlngth_end =  900.7911 #900.8071  #898.274

    wvlength_range = np.linspace(wvlngth_start, wvlngth_end, len(img[1, :]))

    for j in range(len(y_index_range)):
        spectrum_sum = spectrum_sum + img[row_index_setting + y_index_range[j], :]
        bckgrnd_sum_1 = bckgrnd_sum_1 + bckgrnd_img_1[row_index_setting + y_index_range[j], :]
        bckgrnd_sum_2 = bckgrnd_sum_2 + bckgrnd_img_2[row_index_setting + y_index_range[j], :]
        bckgrnd_sum_3 = bckgrnd_sum_3 + bckgrnd_img_3[row_index_setting + y_index_range[j], :]

    bckgrnd_avg = (bckgrnd_sum_1 + bckgrnd_sum_2 + bckgrnd_sum_3)/3
    spectrum_sum = spectrum_sum - bckgrnd_avg 
    spectrum_sum = spectrum_sum/2

    



main_folder_name = '19_08_2024_Test'
bckgrnd_filename = 'background_1200mm_grating_19082024'
bckgrnd_filename_A = 'background_1200mm_grating_19082024_ver_A'   # Change in Set_A at [95, 2]
bckgrnd_filename_B = 'background_1200mm_grating_19082024_ver_B'   # Change in Set_B at [41, 6]

# Parameters for different gratings
row_index_1800mm_grating = 1744
row_index_1200mm_grating = 1788
row_index_setting = row_index_1200mm_grating

N_column = 5
N_row = 200

column_index_list = [1, 2, 3, 4, 5, 6]
row_index_list = np.arange(1, 200, 1)


for column_index in column_index_list:
    
    if column_index == 1: 
        folder_name = main_folder_name + '/Set_A/Spectrometer_ASI'
        bckgrnd_filename = folder_name + '/'+ bckgrnd_filename_A

    elif column_index > 1:
        folder_name = main_folder_name + '/Set_B/Spectrometer_ASI'
        bckgrnd_filename = folder_name + '/'+ bckgrnd_filename_A
    
    for row_index in row_index_list:

        camera_img_filename = folder_name+'QD_Spec_Plot_['+str(row_index)+' '+str(column_index)+']'
        read_spectrum(folder_name, camera_img_filename, bckgrnd_filename)


