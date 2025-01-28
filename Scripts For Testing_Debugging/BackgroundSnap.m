
addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
pause(2)
% initiziling all the needed folders 
%py.qd_data_folder_creation.create_qd_data_directories()

            % Spectroscope plot saving set environment variables 
            terminate(pyenv); % Clear Python 
            setenv('TCL_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tcl8.6")
            setenv('TK_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tk8.6")

            % Inserting proper pathways needed 
            pyLibraryFolder_Scripts = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools"; 
            insert(py.sys.path, int64(0), pyLibraryFolder_Scripts)
            pause(1)

asi_initialization_return = py.asi_func.init_camera();
pause(2)
%py.asi_func.snap_background(asi_initialization_return,0)
%pause(2)
QD_counter = [78 6]; 
pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d]',QD_counter); 

py.asi_func.snap_image(asi_initialization_return,pyrun_file_text_Spec)