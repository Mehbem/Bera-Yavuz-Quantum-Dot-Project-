import pyvisa

rm = pyvisa.ResourceManager()

# List all connected devices

print(rm.list_resources())

# List only USB devices

usb_devices = [device for device in rm.list_resources() if device.startswith('USB')]

print(usb_devices) 