import os 
import time

# Function for creating directories
def create_directory(qd_data_directory_main, directory_name):

	qd_data_directory_path = os.path.join(qd_data_directory_main, directory_name)
	access_rights = 0o755

	if os.path.exists(qd_data_directory_path):

		print("Directory already exists, skipping directory creation")

	if not os.path.exists(qd_data_directory_path):
	
		try:
			os.makedirs(qd_data_directory_path, access_rights)

		except OSError:
			print("************")
			print("Couldn't create directory", qd_data_directory_path)
	
		else:
			print("************") 
			print("Created directory:", qd_data_directory_path) 

#############################################################
# Define today's date string in the format: 'DD_MM_YYYY_Test'
def date_string():

	time_now = time.localtime()

	year = time_now[0]
	month = time_now[1]
	day = time_now[2]
	hour = time_now[3]
	mins = time_now[4]
	sec = time_now[5]

	if (month < 10 and day < 10):
		date = "0"+str(day)+"_0"+str(month)+"_"+str(year)

	elif (month < 10 and day == 10):
		date = str(day)+"_0"+str(month)+"_"+str(year)

	elif (month < 10 and day > 10):   
		date = str(day)+"_0"+str(month)+"_"+str(year)
		
	elif (month > 9 and day < 10):
		date = "0"+str(day)+"_"+str(month)+"_"+str(year)
		
	else: 
		date = str(day)+"_"+str(month)+"_"+str(year)

	return date

#############################################
# Create folders for today's experimental data 
def create_qd_data_directories():

	date = date_string()

	qd_data_directory_main = r"c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data"
	qd_data_directory_testing = r"C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts For Testing_Debugging\Testing Images"
	qd_data_directory_today = date+"_Test/"
	qd_data_directory_today_A = date+"_Test/Position_uEYE/"
	qd_data_directory_today_B = date+"_Test/Spectrometer_ASI/"
	qd_data_directory_today_C = date+"_Test/Spectrometer_Plots/"
	qd_data_director_text_files = date+"_Test/QD_Text_File_Data/"
	qd_data_test_images = date+"_Test/"
	qd_data_directory_Emission_List = date+"_Test/"+"Emission_List"
	qd_data_directory_Plots_FSS = date+"_Test/"+"FSS_Measurements/Spectrometer_Plots_FSS/"

	create_directory(qd_data_directory_main, qd_data_directory_today)
	create_directory(qd_data_directory_main, qd_data_directory_today_A)
	create_directory(qd_data_directory_main, qd_data_directory_today_B)
	create_directory(qd_data_directory_main, qd_data_directory_today_C)
	create_directory(qd_data_directory_testing,qd_data_test_images)
	create_directory(qd_data_directory_main,qd_data_directory_Emission_List)
	create_directory(qd_data_directory_main,qd_data_directory_Plots_FSS)
	create_directory(qd_data_directory_main,qd_data_director_text_files)






