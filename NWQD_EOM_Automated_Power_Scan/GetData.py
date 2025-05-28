import MessangeBox
from ctypes import cdll,c_long, c_ulong, c_uint32,byref,create_string_buffer,c_bool,c_char_p,c_int,c_int16,c_double, sizeof, c_void_p
from TLPM import TLPM

resourceName = None

def initializePm16():
    tlPM = TLPM()
    deviceCount = c_uint32()
    try:
        tlPM.findRsrc(byref(deviceCount))
        print('Device count', deviceCount)
    except Exception:
        MessangeBox.showerror("Error","Please connect PM16 to computer")
        return


    global resourceName
    resourceName = [create_string_buffer(1024), create_string_buffer(1024)]
        
    tlPM.getRsrcName(c_int(0), resourceName[0])
    #rint(c_int(0))
    
    # tlPM.getRsrcName(c_int(1), resourceName[1])
    # print(c_int(1))
    
    tlPM.close()
    pass

def startGetData():
    i = 0 
    
    try:
        global resourceName
        tlPM = TLPM()
        tlPM.open(resourceName[i], c_bool(True), c_bool(True))

        power = c_double()
        tlPM.measPower(byref(power))
        intensity = power.value
        tlPM.close()
        return intensity
    except Exception:
        MessangeBox.showerror("","Error while getting data from PM16")
        return startGetData()

