import numpy as np
import sys 
import qd_data_folder_creation 

print("Date", qd_data_folder_creation.date_string())

print("********************************************")
print("Starting: Creating Folders for Today")

qd_data_folder_creation.create_qd_data_directories() 


