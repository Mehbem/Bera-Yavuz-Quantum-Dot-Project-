import time
import elliptec

addresses = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',\
             'A', 'B', 'C', 'D', 'E', 'F']

def Initialize_ELL14(port1=None, port2=None):
    if(port1 != None):
        controller1 = elliptec.Controller(port1, debug=False)
        #input('Connect characterizing QWP and press Enter')
        for i in addresses:
            try:
                characterize_QWP = elliptec.Rotator(controller1, address=i)
                break
            except:
                continue
        characterize_QWP.change_address('F')

        input('Connect characterizing HWP and press Enter')
        #controller1.close()
        #controller1 = elliptec.Controller(port1)
        for i in addresses:
            try:
                characterize_HWP = elliptec.Rotator(controller1, address=i)
                break
            except:
                continue
        sw = input('Swap?')
        if(sw == 'y'):
            temp = characterize_QWP
            characterize_QWP = characterize_HWP
            characterize_HWP = temp
        characterize_HWP.change_address('1')
        characterize_QWP.change_address('0')
        
    else:
        characterize_QWP = None
        characterize_HWP = None


    if(port2 !=None):
        controller0 = elliptec.Controller(port2, debug=False)
        #input('Connect rotating QWP and press Enter')
        for i in addresses:
            try:
                rotate_QWP = elliptec.Rotator(controller0, address=i)
                break
            except:
                continue
        rotate_QWP.change_address('F')

        #input('Connect rotating HWP and press Enter')
        #controller0.close()
        #controller0 = elliptec.Controller(port2)
        for i in addresses:
            try:
                rotate_HWP = elliptec.Rotator(controller0, address=i)
                break
            except:
                continue
        sw = input('Swap?')
        if(sw == 'y'):
            temp = rotate_QWP
            rotate_QWP = rotate_HWP
            rotate_HWP = temp
        rotate_HWP.change_address('1')
        rotate_QWP.change_address('0')
    else:
        rotate_QWP = None
        rotate_HWP = None
    
    return characterize_QWP, characterize_HWP, rotate_QWP, rotate_HWP
