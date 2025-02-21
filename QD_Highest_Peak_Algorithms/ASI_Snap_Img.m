function [Emission_Reading_Img,pks,plot_img] = ASI_Snap_Img(obj,vid_ASI,src_ASI,ImgType,SaveImg,Spectrometer_Gratting,QD_ID,FSS)
        % ASI_Snap_Img - Captures and processes an image from the ASI Spectrometer.
        %
        % Syntax:
        %   [Emission_Reading_Img] = ASI_Snap_Img(obj, vid_ASI, src_ASI, ImgType, SaveImg, Spectrometer_Gratting, QD_ID)
        %
        % Description:
        %   This function captures an image from the ASI spectrometer system, processes the image based on the 
        %   selected image type ("Background" or "Spectrometer"), and applies necessary corrections such as 
        %   background subtraction and grayscale conversion. If the image type is "Spectrometer," it extracts 
        %   spectral data, detects peak intensities, and saves a wavelength spectrum plot.
        %
        % Inputs:
        %   obj                  - Object instance (not explicitly used in the function but required by class methods)
        %   vid_ASI              - Video input object for the ASI camera
        %   src_ASI              - Source properties of the ASI camera
        %   ImgType              - Type of image to capture: "Background" or "Spectrometer"
        %   SaveImg              - Boolean flag to save the captured image (not explicitly used in the function)
        %   Spectrometer_Gratting - Spectrometer grating setting (e.g., 1800, 1200, 150)
        %   QD_ID                - Identifier for the Quantum Dot sample being analyzed
        %   FSS                     - Specifies its for fine struct splitting (FSS degrees #)
        %
        % Outputs:
        %   Emission_Reading_Img  - Captured and processed emission image (grayscale)
            data = load("Spectrometer_Settings.mat");

            date = Fetch_Date(obj); % fetching today's date
            date_test =  date + "_Test"; 
            
            % Main pathway for where the background photos end up 
            pathway_main = "c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\" + date_test +"\Spectrometer_ASI"; 
            filename_background = sprintf("background_%dmm_grating_%s",Spectrometer_Gratting,date);




            switch ImgType
                case "Background"
                    
                    % Define background image names and storage path
                    num_background_images = 3;

                    % Capture and save background images
                    start(vid_ASI)
                    for i = 1:num_background_images
                        filename_background_saved = sprintf("%s_%d.png", filename_background, i);
                        full_pathway = fullfile(pathway_main, filename_background_saved);
                        background_img(:,:,:,i) = getsnapshot(vid_ASI); 
                        imwrite(background_img(:,:,:,i), full_pathway); % Save image
                        fprintf("Saved background image: %s\n", filename_background_saved);
                    end
                    stop(vid_ASI)
                    

                    % get average of background and store it 
                    background_img = mean(background_img, 4); 
                    save("Spectrometer_Settings.mat","background_img","-append")


                case "Spectrometer"

                    %Parameters for different gratings
                    if Spectrometer_Gratting == 1800

                        % Define Wavelength Range
                        wvlngth_start = data.min_wvlen_1800mm;
                        wvlngth_end = data.max_wvlen_1800mm; 
                    elseif Spectrometer_Gratting == 1200

                    elseif Spectrometer_Gratting == 150

                    end


                    
                    % snapping a photo and gray scaling it 
                    Emission_Reading_Img = [];
                    num_frames_to_use = 3; 
                    start(vid_ASI)
                    for i = 1:num_frames_to_use
                        Emission_Reading_Img(:,:,:,i) = getsnapshot(vid_ASI);  % Collect frames (ensure this is the correct function)
                    end
                    Emission_Reading_Img = round(mean(Emission_Reading_Img, 4));
                    if size(Emission_Reading_Img,3) == 3 % checks if image is rgb
                        Emission_Reading_Img = rgb2gray(Emission_Reading_Img);
                    end
                    stop(vid_ASI)
                    % grab background average
                    background_Img = data.background_img; 
      
                    
                    % Auto-size vertical window
                    [height,width] = size(Emission_Reading_Img); 
                    central_row = height/2; 
                    window_size = round(height*0.4); % 40% below and above the central_row 
                    valid_rows = max(1, central_row - window_size):min(height, central_row + window_size);

                    
                    % Sum spectrum within window
                    spectrum_sum = sum(Emission_Reading_Img(valid_rows, :), 1);
                    background_Img = sum(background_Img(valid_rows, :), 1);
                    
                    spectrum_sum = spectrum_sum - background_Img ;
                    
    
                    % Init wavelength
                    wvlength = linspace(wvlngth_start, wvlngth_end, size(Emission_Reading_Img, 2));
                   

                    % Define the directory path for saving the plot
                    ASI_plots_directory = strcat(date_test, '\Spectrometer_Plots');
                    Quantum_Dot_Named_File = sprintf("[%d %d]_%dmm_gratting",QD_ID,Spectrometer_Gratting);
                    qd_data_ASI_plots_directory = fullfile("c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data", ASI_plots_directory,Quantum_Dot_Named_File);
                    
                    % Checking if FSS is involved 
                    if contains(FSS,"FSS")
                        FSS_Text = strsplit(FSS); 
                        FSS_Angle = FSS_Text{2}; 
                        latestFolder = find_latestFolder(obj,QD_ID);
                        Quantum_Dot_Named_File = sprintf("[%d %d]_%dmm_gratting_FSS_%s_degrees",QD_ID,Spectrometer_Gratting,FSS_Angle);
                        qd_data_ASI_plots_directory = fullfile(latestFolder,Quantum_Dot_Named_File);
                    end

                   % Evaluating the boolean to see if user wants to save raw image
                    if SaveImg == "Yes"
                        if contains(FSS,"FSS")
                            mkdir(fullfile(latestFolder,"Raw_Imgs"))
                            RawImg_Fullpathway = fullfile(latestFolder,"Raw_Imgs",strcat(Quantum_Dot_Named_File,".png"));
                        else
                            specific_filename_img = sprintf("[%d %d]_%dmm_gratting_ASI_RawImg.png",QD_ID,Spectrometer_Gratting);
                            RawImg_Fullpathway = fullfile(pathway_main,specific_filename_img);
                        end
                        imwrite(Emission_Reading_Img,RawImg_Fullpathway)
                    end
                    
                    
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

                    % text file creation 
                    data_matrix = [wvlength;spectrum_sum];
                    specific_filename = sprintf("\\QD_Text_File_Data\\[%d %d]_%dmm_gratting.txt",QD_ID,Spectrometer_Gratting); 
                    text_file_pathway = fullfile("c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\" + date_test,specific_filename);

                    % Checking if FSS is involved 
                    if contains(FSS,"FSS")                      
                        Quantum_Dot_Named_File = sprintf("[%d %d]_%dmm_gratting_FSS_%s_degrees.txt",QD_ID,Spectrometer_Gratting,FSS_Angle);
                        text_file_pathway = fullfile(latestFolder,Quantum_Dot_Named_File);
                    end

                    % opening the text file 
                    file = fopen(text_file_pathway,'w');
                    
                    % Print the settings
                    fprintf(file,'-----Settings-----\n');
                    fprintf(file,'Spectrometer Grating: %dmm \n',Spectrometer_Gratting);
                    fprintf(file,'Exposure: \t %.4E seconds \n', src_ASI.Exposure);
                    fprintf(file,'Gain: \t %.1f \n',src_ASI.Gain);
                    fprintf(file,'Resolution: \t %.0f x %.0f \n', width,height);
                    fprintf(file,'Captured time: %s\n', datestr(now,'HH:MM:SS mmm-dd-yyyy'));
                    % Print the data
                    fprintf(file,'-----Data-----\n');
                    fprintf(file,'%s \t %s\n', 'Wavelen. (nm)','Background Column Sum');
                    fprintf(file,'%.4f \t %u\n',data_matrix);


            end