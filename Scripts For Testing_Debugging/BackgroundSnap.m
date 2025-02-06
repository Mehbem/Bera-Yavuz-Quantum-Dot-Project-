
addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB");

date = MyFuncs.Fetch_Date(); % fetching today's date
date_test =  date + "_Test"; 

num_background_images = 3;
background_images = cell(1, num_background_images); % Preallocate cell array

 % Assigning the default settings to ASI
            ASI_Settings.desired_exposure_time = 2;
            ASI_Settings.Gain = 300;
            ASI_Settings.Gamma = 50; 
            ASI_Settings.Brightness = 50; 

            % Assigning the default settings to UI
            UI_Settings = ""; % no settings as of yet due to only require basic snapping function 

      info = imaqhwinfo;
        if isempty(info.InstalledAdaptors)
            error('No cameras found. Connect a camera and try again.');
        end
        adaptor = info.InstalledAdaptors{1};
        camInfo = imaqhwinfo(adaptor);
        numCameras = numel(camInfo.DeviceIDs); 

        % Check if at least two cameras are connected (making sure the UI and ASI
        % camera are there) 
        if numCameras < 2
            error('Two cameras are required but only %d detected.', numCameras);
        end
         
        % Identify each camera
        for i = 1:numCameras
            deviceID = camInfo.DeviceIDs{i}; % Get device ID
            deviceName = camInfo.DeviceInfo(i).DeviceName; % Get device name
            fprintf('Camera %d: ID = %d, Name = %s\n', i, deviceID, deviceName);
            if contains(deviceName,"ASI") %  checks to see if the detected device is the ASI camera
                ASI_Device_ID = i;
                camInfo.DeviceInfo(i).DefaultFormat = 'RGB8_6248x4176';%'RGB8_1280x960'; 
                fprintf("found ASI (spectrometer camera)\n")
            elseif contains(deviceName,"UI") %  checks to see if the detected device is the ASI camera
                UI_Device_ID = i;
                fprintf("found UI 148x Camera (nanowire camera)\n")
            else % any other camera would have to be an unknown device
                fprintf("unknown device detected\n")
            end
        end
        
        % establishing connection and parameters of ASI Device 
        vid_ASI = videoinput(adaptor, ASI_Device_ID,); % function can take a third input to specify formatting (ASK Sreesh) 
        src_ASI = getselectedsource(vid_ASI);
        %all_props_ASI = propinfo(vid_ASI); shows all settings user can change 
        
        % establishing connection and parameters of UI Device 
        vid_UI = videoinput(adaptor, UI_Device_ID);
        src_UI = getselectedsource(vid_UI);
        % all_props_UI = propinfo(vid_UI); % shows all settings user can change 
        

           all_resolutions = camInfo.DeviceInfo(2).SupportedFormats;
           camInfo.DeviceInfo(2).DefaultFormat = 'RGB8_6248x4176'; 
           
            % Vertically Flipping the UI camera to the proper orientation   
            src_UI.VerticalFlip = 'on'; 

           src_UI.ExposureMode = "manual"; 
            src_UI.GainMode = "manual";
            src_UI.ContrastMode = "manual";

           vid_UI.ReturnedColorSpace = 'grayscale'; % rgb, grayscale, bayer
           vid_UI.VideoFormat = 'RGB8_6248x4176'
% snapping a photo and gray scaling it 
                    Emission_Reading_Img = getsnapshot(vid_ASI);
                    if size(Emission_Reading_Img,3) == 3 % checks if image is rgb
                        Emission_Reading_Img = rgb2gray(Emission_Reading_Img);
                    end

                        QD_ID = [1 1];
                        Spectrometer_Gratting = 1800;
                       pathway_main = "c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\" + date_test +"\Spectrometer_ASI"; 
                        specific_filename_img = sprintf("[%d %d]_%dmm_gratting_ASI_RawImg.png",QD_ID,Spectrometer_Gratting);
                        RawImg_Fullpathway = fullfile(pathway_main,specific_filename_img);
                        imwrite(Emission_Reading_Img,RawImg_Fullpathway)
 
                    
                    %Highest-pixel based threshold:
                    auto_thresh = max(Emission_Reading_Img,[],"all")*0.3; 
                    
                    % Apply threshold
                    filtered_img = Emission_Reading_Img;
                    
                    filtered_img(Emission_Reading_Img <= auto_thresh) = 0;
                                filename_background = sprintf("background_%dmm_grating_%s",Spectrometer_Gratting,date);
                    imshow(filtered_img)
% Read Background Images in grayscale
                    for i = 1:num_background_images
                    filename_background_saved = sprintf("%s_%d.png", filename_background, i);
                    full_pathway = fullfile(pathway_main, filename_background_saved);
                    img = imread(full_pathway);
                    
                        if size(img, 3) == 3 % Convert RGB images to grayscale
                            img = rgb2gray(img);
                        end
                    
                    background_images{i} = img; % Store in cell array
                    
                    end
                    
                    % assigning backgrounds to proper variables 
                    bckgrnd_img_1 = background_images{1};
                    bckgrnd_img_2 = background_images{2};
                    bckgrnd_img_3 = background_images{3};


                   
                    % Find central row with maximum intensity
                    central_row = size(Emission_Reading_Img,1)/2;
                    
                    % Auto-size vertical window
                    [height,width] = size(img); 
                    window_size = round(height*0.3); % 15% below and above the central_row 
                    valid_rows = max(1, central_row - window_size):min(height, central_row + window_size);

                    
                    % Sum spectrum with automated window
                    spectrum_sum = sum(double(Emission_Reading_Img(valid_rows, :)), 1);
                    bckgrnd_sum_1 = sum(double(bckgrnd_img_1(valid_rows, :)), 1);
                    bckgrnd_sum_2 = sum(double(bckgrnd_img_2(valid_rows, :)), 1);
                    bckgrnd_sum_3 = sum(double(bckgrnd_img_3(valid_rows, :)), 1);
                    
                    % averaging out and subtracting background
                    bckgrnd_avg = (bckgrnd_sum_1 + bckgrnd_sum_2 + bckgrnd_sum_3)/3;
                    % saving background_avg into a mat file to use later 
                    save("Spectrometer_Settings.mat","bckgrnd_avg","-append"); 
                    
                    spectrum_sum = spectrum_sum - bckgrnd_avg ;
                    spectrum_sum = spectrum_sum/2;
    
                    % Init wavelength
                    wvlngth_start = 889.1;
                    wvlngth_end = 896;
                    wvlength = linspace(wvlngth_start, wvlngth_end, size(filtered_img, 2));
                   

                    % Define the directory path for saving the plot
                    ASI_plots_directory = strcat(date_test, '\Spectrometer_Plots');
                    Quantum_Dot_Named_File = sprintf("[%d %d]_%dmm_gratting",QD_ID,Spectrometer_Gratting);
                    qd_data_ASI_plots_directory = fullfile("c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data", ASI_plots_directory,Quantum_Dot_Named_File);
                   
                    % Create a new figure
                    figure;
                    
                    % Make the figure invisible
                    set(gcf, 'Visible', 'off');
                    
                    
                    % Plot the spectrum data
                    plot(wvlength, spectrum_sum);
                    
                    % Set title and labels
                    title_font = sprintf("QD Spectrum Plot: [%d %d]",QD_ID);
                    title(title_font);
                    xlabel('Wavelength [nm]');
                    ylabel('Arb. Counts');

                    % plot image
                    plot_img = gcf; 
                    
                    % Adjust the figure size
                    set(gcf, 'Position', [100, 100, 1200, 800]); % [left bottom width height]
        
                    % finding the main peaks 
                    [pks,locs,~,~] = findpeaks(spectrum_sum,wvlength,'SortStr','descend','NPeaks',3,'MinPeakDistance',0.5);

                    % Create text strings for top 3 peaks
                    peak_text = cell(3,1);
                    for n = 1:min(3,length(pks))
                        peak_text{n} = sprintf('Peak %d: %.3f nm Abs. Counts: %.2f', n, locs(n),pks(n));
                    end
                    
                    % Add text box in top-right corner
                    text(0.95, 0.95, peak_text,...
                        'Units', 'normalized',...
                        'HorizontalAlignment', 'right',...
                        'VerticalAlignment', 'top',...
                        'BackgroundColor', [1 1 1 0.7],... % Semi-transparent white
                        'EdgeColor', 'k',...
                        'FontSize', 10,...
                        'Margin', 3);
                    
                    % Set title and labels
                    title_font = sprintf("QD Spectrum Plot: [%d %d]",QD_ID);
                    title(title_font);
                    xlabel('Wavelength [nm]');
                    ylabel('Arb. Counts');


                    % Save the plot as a .png file
                    saveas(gcf, strcat(qd_data_ASI_plots_directory, '_wvl_graph.png'));

                    imaqreset