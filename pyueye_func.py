from pyueye import ueye
import ctypes
import numpy as np
import cv2
import sys
import os
import time 
import qd_data_folder_creation

# Global variable for camera initialization parameters
pyueye_initialization_return = None


def init_camera():
    
    #---------------------------------------------------------------------------------------------------------------------------------------
    #Variables
    hCam = ueye.HIDS(0)             #0: first available camera;  1-254: The camera with the specified camera ID
    sInfo = ueye.SENSORINFO()
    cInfo = ueye.CAMINFO()
    pcImageMemory = ueye.c_mem_p()
    MemID = ueye.int()
    rectAOI = ueye.IS_RECT()
    pitch = ueye.INT()
    nBitsPerPixel = ueye.INT(24)    #24: bits per pixel for color mode; take 8 bits per pixel for monochrome
    channels = 3                    #3: channels for color mode(RGB); take 1 channel for monochrome
    m_nColorMode = ueye.INT()		# Y8/RGB16/RGB24/REG32
    bytes_per_pixel = int(nBitsPerPixel / 8)
    #---------------------------------------------------------------------------------------------------------------------------------------
    print("START")

    # Starts the driver and establishes the connection to the camera
    nRet = ueye.is_InitCamera(hCam, None)
    if nRet != ueye.IS_SUCCESS:
        print("is_InitCamera ERROR")

    # Reads out the data hard-coded in the non-volatile camera memory and writes it to the data structure that cInfo points to
    nRet = ueye.is_GetCameraInfo(hCam, cInfo)
    if nRet != ueye.IS_SUCCESS:
        print("is_GetCameraInfo ERROR")

    # You can query additional information about the sensor type used in the camera
    nRet = ueye.is_GetSensorInfo(hCam, sInfo)
    if nRet != ueye.IS_SUCCESS:
        print("is_GetSensorInfo ERROR")

    nRet = ueye.is_ResetToDefault( hCam)
    if nRet != ueye.IS_SUCCESS:
        print("is_ResetToDefault ERROR")

    # Set display mode to DIB
    nRet = ueye.is_SetDisplayMode(hCam, ueye.IS_SET_DM_DIB)

    # Set the right color mode
    if int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_BAYER:
        # setup the color depth to the current windows setting
        ueye.is_GetColorDepth(hCam, nBitsPerPixel, m_nColorMode)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("IS_COLORMODE_BAYER: ", )
        print("\tm_nColorMode: \t\t", m_nColorMode)
        print("\tnBitsPerPixel: \t\t", nBitsPerPixel)
        print("\tbytes_per_pixel: \t\t", bytes_per_pixel)
        print()

    elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_CBYCRY:
        # for color camera models use RGB32 mode
        m_nColorMode = ueye.IS_CM_BGRA8_PACKED
        nBitsPerPixel = ueye.INT(32)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("IS_COLORMODE_CBYCRY: ", )
        print("\tm_nColorMode: \t\t", m_nColorMode)
        print("\tnBitsPerPixel: \t\t", nBitsPerPixel)
        print("\tbytes_per_pixel: \t\t", bytes_per_pixel)
        print()

    elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_MONOCHROME:
        # for color camera models use RGB32 mode
        m_nColorMode = ueye.IS_CM_MONO8
        nBitsPerPixel = ueye.INT(8)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("IS_COLORMODE_MONOCHROME: ", )
        print("\tm_nColorMode: \t\t", m_nColorMode)
        print("\tnBitsPerPixel: \t\t", nBitsPerPixel)
        print("\tbytes_per_pixel: \t\t", bytes_per_pixel)
        print()

    else:
        # for monochrome camera models use Y8 mode
        m_nColorMode = ueye.IS_CM_MONO8
        nBitsPerPixel = ueye.INT(8)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("else")

    # Can be used to set the size and position of an "area of interest"(AOI) within an image
    nRet = ueye.is_AOI(hCam, ueye.IS_AOI_IMAGE_GET_AOI, rectAOI, ueye.sizeof(rectAOI))
    if nRet != ueye.IS_SUCCESS:
        print("is_AOI ERROR")

    width = rectAOI.s32Width
    height = rectAOI.s32Height

    # Prints out some information about the camera and the sensor
    print("Camera model:\t\t", sInfo.strSensorName.decode('utf-8'))
    print("Camera serial no.:\t", cInfo.SerNo.decode('utf-8'))
    print("Maximum image width:\t", width)
    print("Maximum image height:\t", height)
    print()

    #---------------------------------------------------------------------------------------------------------------------------------------

    # Allocates an image memory for an image having its dimensions defined by width and height and its color depth defined by nBitsPerPixel
    nRet = ueye.is_AllocImageMem(hCam, width, height, nBitsPerPixel, pcImageMemory, MemID)
    if nRet != ueye.IS_SUCCESS:
        print("is_AllocImageMem ERROR")
    else:
        # Makes the specified image memory the active memory
        nRet = ueye.is_SetImageMem(hCam, pcImageMemory, MemID)
        if nRet != ueye.IS_SUCCESS:
            print("is_SetImageMem ERROR")
        else:
            # Set the desired color mode
            nRet = ueye.is_SetColorMode(hCam, m_nColorMode)


    # Activates the camera's live video mode (free run mode)
    nRet = ueye.is_CaptureVideo(hCam, ueye.IS_DONT_WAIT)
    if nRet != ueye.IS_SUCCESS:
        print("is_CaptureVideo ERROR")

    # Enables the queue mode for existing image memory sequences
    nRet = ueye.is_InquireImageMem(hCam, pcImageMemory, MemID, width, height, nBitsPerPixel, pitch)
    if nRet != ueye.IS_SUCCESS:
        print("is_InquireImageMem ERROR")
    else:
        print("Press q to leave the programm")

    pyueye_initialization_return = hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet

    #time.sleep(2)

    return pyueye_initialization_return


def snap_image(pyueye_initialization_return, filename): 
    
    # Fetch today's folder string 
    date_string = str(qd_data_folder_creation.date_string())
    date_string_test = date_string+"_Test"

    # Enter directory where the images are being saved
    qd_data_directory = r"c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data"
    uEYE_raw_directory = str(date_string_test)+"\Position_uEYE"
    qd_data_uEYE_raw_directory = os.path.join(qd_data_directory, uEYE_raw_directory)
    os.chdir(qd_data_uEYE_raw_directory) 

    hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet = pyueye_initialization_return

    #time.sleep(2)
    #print('filename', filename)
    #filename = "test-2.jpg"

    # In order to display the image in an OpenCV window we need to...
    # ...extract the data of our image memory
    array = ueye.get_data(pcImageMemory, width, height, nBitsPerPixel, pitch, copy=False)

    # bytes_per_pixel = int(nBitsPerPixel / 8)

    # ...reshape it in an numpy array...
    frame = np.reshape(array,(height.value, width.value, bytes_per_pixel))

    # ...resize the image by a half
    # frame = cv2.resize(frame,(0,0),fx=0.5, fy=0.5)

    #---------------------------------------------------------------------------------------------------------------------------------------
    #...and finally display it
    #cv2.imshow("SimpleLive_Python_uEye_OpenCV", frame)
    cv2.imwrite(filename, frame) 
    #---------------------------------------------------------------------------------------------------------------------------------------

    cv2.destroyAllWindows()
    #print("A single snap was taken")

def snap_image_test(pyueye_initialization_return, filename): 
    
    # Fetch today's folder string 
    date_string = str(qd_data_folder_creation.date_string())
    date_string_test = date_string+"_Test"

    # Enter directory where the images are being saved
    qd_data_directory = r"C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts For Testing_Debugging\Testing Images"
    uEYE_raw_directory = str(date_string_test)
    qd_data_uEYE_raw_directory = os.path.join(qd_data_directory, uEYE_raw_directory)
    os.chdir(qd_data_uEYE_raw_directory) 

    hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet = pyueye_initialization_return

    #time.sleep(2)
    #print('filename', filename)
    #filename = "test-2.jpg"

    # In order to display the image in an OpenCV window we need to...
    # ...extract the data of our image memory
    array = ueye.get_data(pcImageMemory, width, height, nBitsPerPixel, pitch, copy=False)

    # bytes_per_pixel = int(nBitsPerPixel / 8)

    # ...reshape it in an numpy array...
    frame = np.reshape(array,(height.value, width.value, bytes_per_pixel))

    # ...resize the image by a half
    # frame = cv2.resize(frame,(0,0),fx=0.5, fy=0.5)

    #---------------------------------------------------------------------------------------------------------------------------------------
    #...and finally display it
    #cv2.imshow("SimpleLive_Python_uEye_OpenCV", frame)
    cv2.imwrite(filename, frame) 
    #---------------------------------------------------------------------------------------------------------------------------------------

    cv2.destroyAllWindows()
    #print("A single snap was taken")


def led_snap_image(pyueye_initialization_return, filename): 
    
    # Directory for saving pictures
    directory = r"C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images"
    os.chdir(directory) 

    hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet = pyueye_initialization_return

    #time.sleep(2)
    #print('filename', filename)
    #filename = "test-2.jpg"

    # In order to display the image in an OpenCV window we need to...
    # ...extract the data of our image memory
    array = ueye.get_data(pcImageMemory, width, height, nBitsPerPixel, pitch, copy=False)

    # bytes_per_pixel = int(nBitsPerPixel / 8)

    # ...reshape it in an numpy array...
    frame = np.reshape(array,(height.value, width.value, bytes_per_pixel))

    # ...resize the image by a half
    # frame = cv2.resize(frame,(0,0),fx=0.5, fy=0.5)

    #---------------------------------------------------------------------------------------------------------------------------------------
    #...and finally display it
    #cv2.imshow("SimpleLive_Python_uEye_OpenCV", frame)
    cv2.imwrite(filename, frame) 
    #---------------------------------------------------------------------------------------------------------------------------------------

    cv2.destroyAllWindows()

    #print("A single snap was taken")

def live_video(pyueye_initialization_return):


    hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet = pyueye_initialization_return

    # Continuous image display
    while(nRet == ueye.IS_SUCCESS):

        # In order to display the image in an OpenCV window we need to...
        # ...extract the data of our image memory
        array = ueye.get_data(pcImageMemory, width, height, nBitsPerPixel, pitch, copy=False)

        # bytes_per_pixel = int(nBitsPerPixel / 8)

        # ...reshape it in an numpy array...
        frame = np.reshape(array,(height.value, width.value, bytes_per_pixel))

        # ...resize the image by a half
        frame = cv2.resize(frame,(0,0),fx=0.5, fy=0.5)
        
    #---------------------------------------------------------------------------------------------------------------------------------------
        #Include image data processing here

    #---------------------------------------------------------------------------------------------------------------------------------------

        #...and finally display it
        cv2.imshow("SimpleLive_Python_uEye_OpenCV", frame)

        # Press q if you want to end the loop
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break
    #---------------------------------------------------------------------------------------------------------------------------------------

def live_snap_video(pyueye_initialization_return, filename, max_time):

    hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet = pyueye_initialization_return

    start_time = time.time() 

    # Continuous image display
    while(nRet == ueye.IS_SUCCESS):

        # Elapsed time         
        elapsed_time = time.time() - start_time

        # In order to display the image in an OpenCV window we need to...
        # ...extract the data of our image memory
        array = ueye.get_data(pcImageMemory, width, height, nBitsPerPixel, pitch, copy=False)

        # bytes_per_pixel = int(nBitsPerPixel / 8)

        # ...reshape it in an numpy array...
        frame = np.reshape(array,(height.value, width.value, bytes_per_pixel))

        # ...resize the image by a half
        frame = cv2.resize(frame,(0,0),fx=0.5, fy=0.5)

    #---------------------------------------------------------------------------------------------------------------------------------------
        #Include image data processing here
    #---------------------------------------------------------------------------------------------------------------------------------------
        # Loop end conditions
        #...and finally display it
        cv2.imshow("SimpleLive_Python_uEye_OpenCV", frame)

        # Press q if you want to end the loop
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

        if elapsed_time > max_time:
            cv2.imwrite(filename, frame)  
            break 
    #---------------------------------------------------------------------------------------------------------------------------------------


def exit_camera(pyueye_initialization_return):

    hCam, pcImageMemory, width, height, nBitsPerPixel, pitch, bytes_per_pixel, MemID, nRet = pyueye_initialization_return
        
    # Releases an image memory that was allocated using is_AllocImageMem() and removes it from the driver management
    ueye.is_FreeImageMem(hCam, pcImageMemory, MemID)

    # Disables the hCam camera handle and releases the data structures and memory areas taken up by the uEye camera
    ueye.is_ExitCamera(hCam)

    # Destroys the OpenCv windows
    cv2.destroyAllWindows()

    print("Camera has exited")
    print("END")

# t_start = time.time()
# pyueye_initialization_return = init_camera()
# t_end = time.time()

# snap_image(pyueye_initialization_return, 'success_ueye.jpg')
# # print('Time for initializing camera', (t_end - t_start))

# # live_video(pyueye_initialization_return)

# exit_camera(pyueye_initialization_return)
