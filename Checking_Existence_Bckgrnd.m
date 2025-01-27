% Bera Yavuz 
clear;
clc;

           % adding path way for function container and python scripts used
            addpath("C:\Users\Quantum Dot\Desktop\Bera_Yavuz_GitHub\AttoCube-Project-Stuff","C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
        
            % adding pathway for file with images for finding excitation laser
            addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images")

            %Fetch today's folder string 
            date = py.qd_data_folder_creation.date_string(); 
            date = string(date);

pathway_name_original = sprintf("C:\\Users\\Quantum Dot\\Desktop\\Bera Yavuz - ANC300 Movement and Images\\QD_Data\\%s_Test\\Spectrometer_ASI\\background_1800mm_grating_%s_1.png",date,date);
pathway_name_FSS = sprintf("C:\\Users\\Quantum Dot\Desktop\\Bera Yavuz - ANC300 Movement and Images\\QD_Data\\%s_Test\\FSS_Measurements\\Spectrometer_ASI_FSS",date);


if ~exist(pathway_name_original,'file') | ~exist(pathway_name_FSS,'file')
    fprintf("background photo has not been taken for the day. Please take a photo before proceeding")
    return
end

fprintf("hello")