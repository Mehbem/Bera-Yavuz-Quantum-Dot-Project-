% Snapping Photo and saving it with today's date
            % Filtering Settings 
            scaling = 0.5; 
            sigma_flatfield = 60; 
            Salt_pepper_pixel_factor = [5 5]; 

            [UI_Position_Img] = MyFuncs.UI_Snap_Img(vid_UI,src_UI,"LED_Background_Debug",""); 

            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            figure()
                     
            % Scailing Image for faster processing  
            grayImage = imresize(UI_Position_Img, scaling);
            
            % Store the original dimensions
            [originalHeight, originalWidth, ~] = size(UI_Position_Img);
            
            % Making the shadow across the image more balanced 
            FlatFieldAdjustedImg = imflatfield(grayImage,sigma_flatfield);
            imshow(FlatFieldAdjustedImg)

            % Getting rid of the salt-pepper noise 
            imgFiltered = medfilt2(FlatFieldAdjustedImg,Salt_pepper_pixel_factor);
            imshow(imgFiltered)
            % Binarizing Image for simplicity 
            LED_Spot = imbinarize(imgFiltered,"adaptive","ForegroundPolarity","bright","Sensitivity",0.3);
            %LED_Spot = bwareaopen(LED_Spot,50);
            
            
            % Filtering Based on circularity 
            labelledImage = bwlabel(LED_Spot);
            imshow(labelledImage)
            propertiesImage = regionprops(labelledImage,"Circularity","Centroid","Area"); 
            AllAreas = [propertiesImage.Area]; 
            AllCircularities = [propertiesImage.Circularity];
            AllCentroids = vertcat(propertiesImage.Centroid); 
            ValidRegions = find((AllCircularities > 0.9) & AllAreas > 100);
            LED_Spot_BW = ismember(labelledImage,ValidRegions);

            LED_centroid = AllCentroids(ValidRegions,:);
            imshow(LED_Spot_BW)

            
 
            % Scailing everything back up 
            LED_centroid = LED_centroid ./scaling;

            % Updating the LED coordinate text file with the new coordinate 
            MyFuncs.LED_Cooordinate_Identifer(LED_centroid,"Write","LAB"); 
            LED_Spot_BW = imresize(LED_Spot_BW,[originalHeight, originalWidth]); 
            
            GUI_orange = [0.85,0.33,0.10]; 
            visibility_setting = "off"; 
            % Plotting everything for visual purposes 
            figure("Name","OriginalImage","NumberTitle","off","Color",GUI_orange,"Visible",visibility_setting);
            imshow(UI_Position_Img)
            title("Original Image")
            exportgraphics(gcf, 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images\LED_Process_plot_Step1.png', 'Resolution', 300, 'BackgroundColor', GUI_orange);
            hold off; 
            
            figure("Name","First Filter","NumberTitle","off","Color",GUI_orange,"Visible",visibility_setting);
            imshow(LED_Spot)
            title("Filtered Image")
            hold on;
            plot(AllCentroids(:,1),AllCentroids(:,2),"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o","LineStyle","none");
            exportgraphics(gcf, 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images\LED_Process_plot_Step2.png', 'Resolution', 300, 'BackgroundColor', GUI_orange);
            hold off; 
            
            figure("Name","Final Filter","NumberTitle","off","Color",GUI_orange,"Visible",visibility_setting);
            title_text = sprintf("The coordinates of the excitation laser is [%.3f,%.3f]",LED_centroid(1),LED_centroid(2));
            imshow(LED_Spot_BW)
            title(title_text)
            hold on;
            plot(LED_centroid(1),LED_centroid(2),"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o");
            exportgraphics(gcf, 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images\LED_Process_plot_Final.png', 'Resolution', 300, 'BackgroundColor', GUI_orange);
            hold off; 
            
            