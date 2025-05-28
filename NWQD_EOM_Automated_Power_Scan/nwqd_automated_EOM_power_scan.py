

# import for communicating with ThorLabs and ELL14K
import GetData

# import socket for communication with SIGLENT waveform generator
import socket
import sys
import time

# import numpy and matplotlib
import matplotlib.pyplot as plt
import numpy as np

# import ASI camera
import asi_func

import os
import qd_data_folder_creation

#ooooooooooooooooooooooooooooo#
#ooooooooooooooooooooooooooooo#
# Initialize Thorlabs power meter
GetData.initializePm16()
print('Initialized PM')

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
asi_initialization_return = asi_func.init_camera() 

# Input QD [row, column] indices
qd_row_val = 88
qd_column_val = 1

# Create folders for that QD 
qd_data_power_scan_data_directory = r"c:\Users\Quantum Dot\Desktop\NWQD_Operation\Scripts_&_Debugging_Tools\NWQD_EOM_Automated_Power_Scan\NWQD_Power_Scan_Data"
bckgrnd_img_directory = 'Background_Images'
qd_specific_directory = 'QD ['+str(qd_row_val)+', '+str(qd_column_val)+']'+'_Power_Scan'
qd_specific_directory_full_path = qd_data_power_scan_data_directory + '\QD ['+str(qd_row_val)+', '+str(qd_column_val)+']'+'_Power_Scan'
qd_specific_ASI_directory = 'ASI_Images'
qd_specific_plot_directory = 'Spectrometer_Plots'

qd_data_folder_creation.create_directory(qd_data_power_scan_data_directory, bckgrnd_img_directory)
qd_data_folder_creation.create_directory(qd_data_power_scan_data_directory, qd_specific_directory)
qd_data_folder_creation.create_directory(qd_specific_directory_full_path, qd_specific_ASI_directory)
qd_data_folder_creation.create_directory(qd_specific_directory_full_path, qd_specific_plot_directory)

bckgrnd_image_directory_full_path = qd_data_power_scan_data_directory + '\Background_Images'
qd_specific_ASI_directory_full_path = qd_specific_directory_full_path + '\ASI_Images'
qd_specific_plot_directory_full_path = qd_specific_directory_full_path + '\Spectrometer_Plots'

#******************************************************************#
#******************************************************************#
# Save background images 

# Do you want to take a new set of measurements? 
new_bckgrnd_images = 0 # Yes/No
if new_bckgrnd_images == 1: 
    user_input = input("Block the excitation beam to take a new set of background images (y/n): ")
    if user_input.lower() == "y":
        print("Background snapping starting...")
    
    asi_func.snap_background(asi_initialization_return, bckgrnd_image_directory_full_path)
    
    user_input = input("Remove the block to proceed (y/n): ")
    if user_input.lower() == "n":
        print("Background snapping complete...")


#******************************************************************#
#******************************************************************#
# Interface with Siglent 
SIGLENT_IP = "10.2.0.11" # should match the instrument’s IP address
SIGLENT_PORT = 5024 # the port number of the instrument service
 
count = 0
 
def SocketConnect():
    try:
        #create an AF_INET, STREAM socket (TCP)
        s_SIGLENT = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    except socket.error:
        print ('Failed to create socket.')
        sys.exit();
    try:
        #Connect to remote server
        s_SIGLENT.connect((SIGLENT_IP , SIGLENT_PORT))
    except socket.error:
        print ('Failed to connect to ip ' + SIGLENT_IP)
    return s_SIGLENT
	
def SocketSend(Sock, cmd):
    try :
        #Send cmd string
        Sock.sendall(cmd)
        Sock.sendall(b'\n')
        time.sleep(1)
    except socket.error:
        #Send failed
        print ('Send failed')
        sys.exit()
    #reply = Sock.recv(4096)
    #return reply
 
def SocketClose(Sock):
    #close the socket
    Sock.close()
    time.sleep(1)

#******************************************************************#
#******************************************************************#
s_SIGLENT = SocketConnect()

N_sample = 20
DC_voltage_sweep = np.linspace(0,5, N_sample)

power_value_log = np.zeros(N_sample)

qStr = SocketSend(s_SIGLENT, b'C1:BSWV WVTP,DC') #Set CH1 Wavetype to DC
qStr = SocketSend(s_SIGLENT, b'C1:BSWV OFST,'+str.encode(str(0.5))) #Set CH1 DC Voltage
qStr = SocketSend(s_SIGLENT, b'C1:OUTP LOAD,50') #Switch ON 
qStr = SocketSend(s_SIGLENT, b'C1:OUTP ON') #Switch ON 
    
# Sweep the Siglent for every fixed setting on the MogLabs
# The Siglent sweep triggers data collection on the SNSPD 
for i in range(N_sample):
    print('**************************')
    print('Iteration No. ', i)
    print('Siglent DC Amp. Set')
    qStr = SocketSend(s_SIGLENT, b'C1:BSWV OFST,'+str.encode(str(DC_voltage_sweep[i]))) #Set CH1 DC Voltage
    time.sleep(0.5)
        
    # Fetch power meter value
    os.chdir(qd_data_power_scan_data_directory)
    power_value_log[i] = GetData.startGetData()
    power_value_log[i] = round(power_value_log[i]*1e6, 2) # Converting to muW
    print('Power measured (muW)', power_value_log[i])

    print('Starting/Saving QD Spectrum')
    asi_func.snap_image(asi_initialization_return, 'QD ['+str(qd_row_val)+', '+str(qd_column_val)+']'+'_Power_Level_'+str(power_value_log[i]), bckgrnd_image_directory_full_path, qd_specific_ASI_directory_full_path, qd_specific_plot_directory_full_path)


SocketClose(s_SIGLENT) #Close socket

print('Complete. Exiting program')
sys.exit

