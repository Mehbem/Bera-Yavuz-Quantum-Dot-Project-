classdef functionsContainer
    
    methods
        
        % Mainly Used in ImageAlgorithmScript and XYANC300Axes
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        function [grayImage,imgFiltered,ReversedBWImg,cleanedImg,filteredMask,filterOutlierMask,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = Analyzed_QDBinaryImg_With_Dots(obj,scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,I)
            % Description: takes raw image and identifies Quantum Dots (QD) alongside their centroids via preprocessing and postprocessing 
            % Inputs:
                % Img - the raw image matrix
                % scaling - used for scaling down image for quicker processing 
                % radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area - parameters used for preprocessing and postprocessing ROI 
            % Outputs: 
                % grayImage,imgFiltered,ReversedBWImg,cleanedImg,filteredMask,filterOutlierMask,Contrast_adjusted_Img,Img,MaskedImage - matrices that return different stages of the image processing 
                % centroid,copycentroid - return a m x 2 matrix consisting of the x and y positions of the QD centroids 
                % Mask, Mask_scaled - logical matrices used for creating the ROI within the frame 

            Img = imread(I);
            Img = im2gray(Img);
            grayImage = imresize(Img, scaling);
            sizeImg = size(grayImage);
            sizeImgScaled = size(Img);


            % Get coordinates of the circle.
            angles = linspace(0, 2*pi, 10000);
            x = cos(angles) * radii_big_circle + center_big_circle(1);
            y = sin(angles) * radii_big_circle + center_big_circle(2);
            % Get a mask of the circle
            mask = poly2mask(x, y, sizeImg(1), sizeImg(2));

            % Scaled Up mask
            x_scaled = x ./ scaling;
            y_scaled = y ./ scaling;
            mask_scaled = poly2mask(x_scaled,y_scaled,sizeImgScaled(1),sizeImgScaled(2));


            % Making the shadowing across the image more equal
            FlatFieldAdjustedImg = imflatfield(grayImage,sigma_flatfield);

            % Define the radius for the thin lining
            liningRadius = 10;


            % Dilate the mask to create the outer lining
            se = strel('disk', liningRadius);
            dilatedMask = imdilate(mask, se);

            % Get the thin outer lining by subtracting the original mask from the dilated mask
            thinLining = dilatedMask & ~mask;

            % Apply mask to clear outside image area
            maskedImage = FlatFieldAdjustedImg;
            maskedImage(~mask) = 0;

            % Making the constrast between the dots and the background more clear
            Contrast_adjusted_Img = maskedImage;
            Contrast_adjusted_Img(mask == 1) = imadjust(maskedImage(mask == 1));
            Contrast_adjusted_Img(thinLining) = 255;

            % Getting Rid of Salt-Pepper Noise
            imgFiltered = medfilt2(Contrast_adjusted_Img,Salt_pepper_pixel_factor);

            % Binarizing the image very specifically
            BwImg = imbinarize(imgFiltered, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.19);

            % Reapply mask outside of circle
            RectanglularBWImg = BwImg;
            RectanglularBWImg(~mask) = 1;

            % Reverse the colour of foreground and background
            ReversedBWImg = ~RectanglularBWImg;
            %imshow(ReversedBWImg)

            % Removing any points or pixels that are anomolies
            cleanedImg = bwareaopen(ReversedBWImg,Min_circle_area);
            %imshow(cleanedImg)


            % Taking the Properties of Remaining Image
            [labeledImage, numObjects] = bwlabel(cleanedImg);
            props = regionprops(labeledImage, 'Area', 'Perimeter', 'Eccentricity',"EulerNumber");
            
    

            % Initialize a mask to keep desired objects
            filteredMask = false(size(cleanedImg));

            % Define thresholds for circularity and eccentricity
            circularityThreshold = 0.5; % Adjust as needed
            eccentricityThreshold = 0.83; % Adjust as needed
            
            % for loop for calculating the circularity
            for k = 1 : numObjects
                area = props(k).Area;
                perimeter = props(k).Perimeter;
                eccentricity = props(k).Eccentricity;
                Hollowness = props(k).EulerNumber; 

                % Calculate circularity
                circularity = 4 * pi * area / (perimeter^2);

                % Filter based on circularity and eccentricity
                if circularity > circularityThreshold && eccentricity < eccentricityThreshold && Hollowness == 1
                    filteredMask(labeledImage == k) = true;
                end
            end

            % Further Filtering Using properties
            [labeledImageV2, numObjectsV2] = bwlabel(filteredMask);
            props_centroid = regionprops(filteredMask, 'Centroid','Area'); % Labelling each centroid
            Centroids = vertcat(props_centroid.Centroid); % Putting them into a list for ease of access
            avg_A_detected_spots = sum([props_centroid.Area])/numObjectsV2; 

            % Doing a minimum distance check to filter out points that have shape and
            % intensity properties but don't follow grid pattern
            filterOutlierMask = false(size(filteredMask));
            for j = 1:numObjectsV2
                DistComp_OutlierQD = sqrt((Centroids(:, 1) - Centroids(j,1)).^2 + (Centroids(:, 2) - Centroids(j,2)).^2);
                SortedDistCompQD = sort(DistComp_OutlierQD);
                MinimumDistCompQD = SortedDistCompQD(2);
                if MinimumDistCompQD <= 80
                    filterOutlierMask(labeledImageV2 == j) = false;
                elseif MinimumDistCompQD <= 180
                    filterOutlierMask(labeledImageV2 == j) = true;
                end
                if props_centroid(j).Area < avg_A_detected_spots/2
                   filterOutlierMask(labeledImageV2 == j) = false; 
                end
            end

            % Array of coordinates scaled and scaling back up
            [FinalbwQDImg,~] = bwlabel(filterOutlierMask);
            statsDots = regionprops(FinalbwQDImg, 'Centroid');
            centroid = vertcat(statsDots.Centroid) ./ scaling;
            CopyCentroid = centroid;
            

             % Masked Image that will be Used for backgrounds
            GrayedImage = im2gray(Img);
            MaskedImage = imadjust(GrayedImage);
            MaskedImage(~mask_scaled) = 0;
        end
        
        function [allNextPts,allPerpPts,p1,p2,p3,p4,p5,p6] = MainAxes(obj,Img,centroid,CopyCentroid,radiusQD) %#ok<*INUSD>
            % Description:
                % - finds a red axis and a blue axis that act as the corresponding x and y axis for the grid pattern of the QD
            % Inputs:
                % Img - raw image matrix
                % centroid, copycentroid - centroid coordinates of the Quantum Dots
                % radiusQD - hardcoded radius of QD for plotting 
            % Outputs:
                % AllNextPts - returns QD on the diagonal of image for red axis (x-axis)
                % AllPerpPts - returns QD on the other diagonal of image for blue axis (y-axis)

            % Calculate distances from each point to the top right corner
            distancesX = sqrt((centroid(:, 1) - size(Img,2)).^2);
            distancesY = sqrt((centroid(:, 2) - 0).^2);
            distances = distancesX + distancesY;

            % Creating a table with the indexs of each point alongside the distance from the top corner
            Indices = 1:length(distances);
            Table_with_dist_ind = table(Indices(:),distances(:),'VariableNames',["Indices", "Distance_to_Corner"]);
            Sorted_Table_with_dist_ind = sortrows(Table_with_dist_ind,"Distance_to_Corner");

            % Estimate the maximum number of points (initial size of coords as an estimate)
            maxPoints = size(centroid, 1);

            % Parameter for Accuracy Wanted
            num_of_QD_diagonal = 5;

            i = 0;
            while i <= length(distances)
                i = i + 1;
                % Get the closest point coordinates
                closestPoint = CopyCentroid(Sorted_Table_with_dist_ind{i,1}, :);
                StartingPt = closestPoint;

                % Finding next closest point
                closestPointNew = closestPoint + [radiusQD, -radiusQD];

                % Preallocate array to store all next points
                allNextPts = zeros(maxPoints+1, 2);


                % For Box pattern creating a lower boundary to prevent skewing of Axes line
                DiameterQD = 2*radiusQD;
                LowerBound = closestPointNew - [2.5*DiameterQD, - 4.5*DiameterQD];
                UpperBound = closestPointNew - [4.5*DiameterQD, - 2.5*DiameterQD];

                % Redefining each set variable
                centroid = CopyCentroid;
                numNextPts = 0;

                % Looping through points to find a diagonal of pts to base all lines off of
                % figure("Name","Axes Finding","Color",skyBlue)
                % imshow(MaskedImage)
                % hold on;

                while ~isempty(centroid)
                    % Testing
                    p1 = plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330],"LineStyle","none");
                    p2 = plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2],"LineStyle","none");
                    % Remove points equal to the closest point
                    targetRow = all(centroid == closestPoint, 2);
                    centroid(targetRow, :) = [];

                    % Check if coords is empty after removal
                    if isempty(centroid)
                        break;
                    end

                    % Making a smaller boxed region below the previous point to search for the next point
                    CoordCompUpper = centroid < UpperBound;
                    CoordCompLower = LowerBound < centroid;
                    % Identify corner points for an imaginary scquare around the next point
                    UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                    LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);

                    % Use logical indexing to keep the point that matches
                    NextPt = centroid(UpperPattern_mask & LowerPattern_mask, :);
                    centroid(UpperPattern_mask & LowerPattern_mask, :) = [];

                    % Add the closest point to the list of all next points
                    if ~isempty(NextPt)
                        numNextPts = numNextPts + 1;
                        CopynumNextPts = numNextPts;
                        allNextPts(numNextPts, :) = NextPt;
                        %plot(NextPt(1),NextPt(2),"Marker","o","MarkerSize",5,"MarkerFaceColor",[0 1 0]);
                    end

                    % Use logical indexing to remove all the points before each lowerbound
                    % to prevent the while loop from going on forever
                    RightOfLowerBound = LowerBound < centroid;
                    PointsRemovedRightSide = (RightOfLowerBound(:,1) == 1);
                    centroid(PointsRemovedRightSide,:) = [];

                    % Use logical indexing to remove all the points above each upperbound
                    % to prevent the while loop from going on forever
                    AboveUpperBound = UpperBound > centroid;
                    PointsRemovedAbove = (AboveUpperBound(:,2) == 1);
                    centroid(PointsRemovedAbove,:) = [];

                    % Check if coords is empty after filtering
                    if isempty(centroid)
                        break;
                    end


                    if isempty(NextPt)
                        LowerBound = closestPointNew - [6*DiameterQD, - 8*DiameterQD];
                        UpperBound = closestPointNew - [8*DiameterQD, - 5*DiameterQD];
                        closestPointNew = (UpperBound + LowerBound)./2;
                        %plot(closestPointNew(1),closestPointNew(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.1 0.1 0.9]); % For testing purposes with new images
                    else
                        % Update the closestPoint, closestPointNew, Lowerbound for the next iteration
                        closestPoint = NextPt;
                        closestPointNew = closestPoint + [radiusQD, -radiusQD];
                        LowerBound = closestPointNew - [2.5*DiameterQD, - 4.5*DiameterQD];
                        UpperBound = closestPointNew - [4.5*DiameterQD, - 2.5*DiameterQD];
                    end


                end
                if numNextPts + 1 >= num_of_QD_diagonal
                    break
                elseif (i == length(distances)) && (numNextPts < num_of_QD_diagonal)
                    num_of_QD_diagonal = num_of_QD_diagonal - 1;
                    i = 0;
                elseif num_of_QD_diagonal == 0
                    error("No possible QD were found on the diagonal")
                end

            end

            % Trim the allNextPts array to the actual number of points added
            allNextPts = vertcat(StartingPt,allNextPts(1:CopynumNextPts, :));

            % Find QD across the Perpendicular diagonal
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            PerpCentroid = CopyCentroid;

            PerpStartingPt = allNextPts(round(length(allNextPts)/2),:);
            PerpLowerBoundDown =  PerpStartingPt + [2*DiameterQD,  4*DiameterQD];
            PerpUpperBoundDown = PerpStartingPt + [4*DiameterQD,  2*DiameterQD];

            PerpLowerBoundUp = PerpStartingPt - [4*DiameterQD, 2*DiameterQD];
            PerpUpperBoundUp = PerpStartingPt - [2*DiameterQD, 4*DiameterQD];


            % Preallocate array to store all next points
            allPerpNextPtsDown = zeros(maxPoints, 2);
            allPerpNextPtsUp = zeros(maxPoints, 2);

            % Counter to trim list
            PerpCounterDown = 0;
            PerpCounterUp = 0;

            while PerpUpperBoundDown(2) <= 1920 || PerpLowerBoundDown(2) <= 2000

                % Downwards Search
                p3 = plot(PerpLowerBoundDown(1),PerpLowerBoundDown(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[1, 0, 1],"LineStyle","none");
                p4 = plot(PerpUpperBoundDown(1),PerpUpperBoundDown(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 1, 0],"LineStyle","none");
                % Making a smaller boxed region to search for the next point
                CoordCompUpperDown = PerpCentroid < PerpUpperBoundDown;
                CoordCompLowerDown = PerpLowerBoundDown < PerpCentroid;
                % Identify corner points for an imaginary scquare around the next point
                UpperPattern_maskDown = (CoordCompUpperDown(:, 1) == 1) & (CoordCompUpperDown(:, 2) == 0);
                LowerPattern_maskDown = (CoordCompLowerDown(:, 1) == 1) & (CoordCompLowerDown(:, 2) == 0);
                % Use logical indexing to keep the point that matches
                PerpNextPtDown = PerpCentroid(UpperPattern_maskDown & LowerPattern_maskDown, :);

                % Checking to see if that Point exists
                if isempty(PerpNextPtDown)
                    PerpNextPtDown = (PerpLowerBoundDown + PerpUpperBoundDown)./2;
                    %plot(PerpNextPtDown(1),PerpNextPtDown(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                    PerpLowerBoundDown = PerpNextPtDown + [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundDown = PerpNextPtDown + [4*DiameterQD,  2*DiameterQD];



                else
                    PerpCounterDown = PerpCounterDown + 1;
                    allPerpNextPtsDown(PerpCounterDown,:) = PerpNextPtDown;
                    PerpLowerBoundDown = PerpNextPtDown + [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundDown = PerpNextPtDown + [4*DiameterQD,  2*DiameterQD];
                    %plot(PerpNextPtDown(1),PerpNextPtDown(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                end

            end

            while PerpUpperBoundUp(2) >= -100 || PerpLowerBoundUp(2) >= 0

                % Upwards Search
                p5 = plot(PerpLowerBoundUp(1),PerpLowerBoundUp(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[1, 0, 1],"LineStyle","none");
                p6 = plot(PerpUpperBoundUp(1),PerpUpperBoundUp(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 1, 0],"LineStyle","none");
                % Making a smaller boxed region to search for the next point
                CoordCompUpperUp = PerpCentroid < PerpUpperBoundUp;
                CoordCompLowerUp = PerpLowerBoundUp< PerpCentroid;
                % Identify corner points for an imaginary scquare around the next point
                UpperPattern_maskUp = (CoordCompUpperUp(:, 1) == 1) & (CoordCompUpperUp(:, 2) == 0);
                LowerPattern_maskUp = (CoordCompLowerUp(:, 1) == 1) & (CoordCompLowerUp(:, 2) == 0);
                % Use logical indexing to keep the point that matches
                PerpNextPtUp = PerpCentroid(UpperPattern_maskUp & LowerPattern_maskUp, :);

                % Checking to see if that Point exists
                if isempty(PerpNextPtUp)
                    PerpNextPtUp = (PerpLowerBoundUp + PerpUpperBoundUp)./2;
                    %plot(PerpNextPtUp(1),PerpNextPtUp(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                    PerpLowerBoundUp = PerpNextPtUp - [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundUp = PerpNextPtUp - [4*DiameterQD,  2*DiameterQD];

                else
                    PerpCounterUp = PerpCounterUp + 1;
                    allPerpNextPtsUp(PerpCounterUp,:) = PerpNextPtUp;
                    PerpLowerBoundUp = PerpNextPtUp - [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundUp = PerpNextPtUp - [4*DiameterQD,  2*DiameterQD];
                    %plot(PerpNextPtUp(1),PerpNextPtUp(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                end
            end


            % Trim the allPerpNextPts array to the actual number of points added
            allPerpPts = vertcat(PerpStartingPt,allPerpNextPtsDown(1:PerpCounterDown, :),allPerpNextPtsUp(1:PerpCounterUp,:));
        end

        function [allNextPts,allPerpPts] = ANC300Axes_MainAxes(obj,Img,centroid,CopyCentroid,radiusQD)
             % Description:
                % - finds a red axis and a blue axis that act as the corresponding x and y axis for the grid pattern of the QD (Does not do plotting like MainAxes function)
            % Inputs:
                % Img - raw image matrix
                % centroid, copycentroid - centroid coordinates of the Quantum Dots
                % radiusQD - hardcoded radius of QD for plotting 
            % Outputs:
                % AllNextPts - returns QD on the diagonal of image for red axis (x-axis)
                % AllPerpPts - returns QD on the other diagonal of image for blue axis (y-axis)

            
            % Calculate distances from each point to the top right corner
            distancesX = sqrt((centroid(:, 1) - size(Img,2)).^2);
            distancesY = sqrt((centroid(:, 2) - 0).^2);
            distances = distancesX + distancesY;

            % Creating a table with the indexs of each point alongside the distance from the top corner
            Indices = 1:length(distances);
            Table_with_dist_ind = table(Indices(:),distances(:),'VariableNames',["Indices", "Distance_to_Corner"]);
            Sorted_Table_with_dist_ind = sortrows(Table_with_dist_ind,"Distance_to_Corner");

            % Estimate the maximum number of points (initial size of coords as an estimate)
            maxPoints = size(centroid, 1);

            % Parameter for Accuracy Wanted
            num_of_QD_diagonal = 5;

            i = 0;
            while i <= length(distances)
                i = i + 1;
                % Get the closest point coordinates
                closestPoint = CopyCentroid(Sorted_Table_with_dist_ind{i,1}, :);
                StartingPt = closestPoint;

                % Finding next closest point
                closestPointNew = closestPoint + [radiusQD, -radiusQD];

                % Preallocate array to store all next points
                allNextPts = zeros(maxPoints+1, 2);


                % For Box pattern creating a lower boundary to prevent skewing of Axes line
                DiameterQD = 2*radiusQD;
                LowerBound = closestPointNew - [2.5*DiameterQD, - 4.5*DiameterQD];
                UpperBound = closestPointNew - [4.5*DiameterQD, - 2.5*DiameterQD];

                % Redefining each set variable
                centroid = CopyCentroid;
                numNextPts = 0;

                while ~isempty(centroid)
                    % Remove points equal to the closest point
                    targetRow = all(centroid == closestPoint, 2);
                    centroid(targetRow, :) = [];

                    % Check if coords is empty after removal
                    if isempty(centroid)
                        break;
                    end

                    % Making a smaller boxed region below the previous point to search for the next point
                    CoordCompUpper = centroid < UpperBound;
                    CoordCompLower = LowerBound < centroid;
                    % Identify corner points for an imaginary scquare around the next point
                    UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                    LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);

                    % Use logical indexing to keep the point that matches
                    NextPt = centroid(UpperPattern_mask & LowerPattern_mask, :);
                    centroid(UpperPattern_mask & LowerPattern_mask, :) = [];

                    % Add the closest point to the list of all next points
                    if ~isempty(NextPt)
                        numNextPts = numNextPts + 1;
                        CopynumNextPts = numNextPts;
                        allNextPts(numNextPts, :) = NextPt;

                    end

                    % Use logical indexing to remove all the points before each lowerbound
                    % to prevent the while loop from going on forever
                    RightOfLowerBound = LowerBound < centroid;
                    PointsRemovedRightSide = (RightOfLowerBound(:,1) == 1);
                    centroid(PointsRemovedRightSide,:) = [];

                    % Use logical indexing to remove all the points above each upperbound
                    % to prevent the while loop from going on forever
                    AboveUpperBound = UpperBound > centroid;
                    PointsRemovedAbove = (AboveUpperBound(:,2) == 1);
                    centroid(PointsRemovedAbove,:) = [];

                    % Check if coords is empty after filtering
                    if isempty(centroid)
                        break;
                    end


                    if isempty(NextPt)
                        LowerBound = closestPointNew - [6*DiameterQD, - 8*DiameterQD];
                        UpperBound = closestPointNew - [8*DiameterQD, - 5*DiameterQD];
                        closestPointNew = (UpperBound + LowerBound)./2;

                    else
                        % Update the closestPoint, closestPointNew, Lowerbound for the next iteration
                        closestPoint = NextPt;
                        closestPointNew = closestPoint + [radiusQD, -radiusQD];
                        LowerBound = closestPointNew - [2.5*DiameterQD, - 4.5*DiameterQD];
                        UpperBound = closestPointNew - [4.5*DiameterQD, - 2.5*DiameterQD];
                    end


                end
                if numNextPts + 1 >= num_of_QD_diagonal
                    break
                elseif (i == length(distances)) && (numNextPts < num_of_QD_diagonal)
                    num_of_QD_diagonal = num_of_QD_diagonal - 1;
                    i = 0;
                elseif num_of_QD_diagonal == 0
                    error("No possible QD were found on the diagonal")
                end

            end

            % Trim the allNextPts array to the actual number of points added
            allNextPts = vertcat(StartingPt,allNextPts(1:CopynumNextPts, :));

            % Find QD across the Perpendicular diagonal
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            PerpCentroid = CopyCentroid;

            PerpStartingPt = allNextPts(round(length(allNextPts)/2),:);
            PerpLowerBoundDown =  PerpStartingPt + [2*DiameterQD,  4*DiameterQD];
            PerpUpperBoundDown = PerpStartingPt + [4*DiameterQD,  2*DiameterQD];

            PerpLowerBoundUp = PerpStartingPt - [4*DiameterQD, 2*DiameterQD];
            PerpUpperBoundUp = PerpStartingPt - [2*DiameterQD, 4*DiameterQD];


            % Preallocate array to store all next points
            allPerpNextPtsDown = zeros(maxPoints, 2);
            allPerpNextPtsUp = zeros(maxPoints, 2);

            % Counter to trim list
            PerpCounterDown = 0;
            PerpCounterUp = 0;

            while PerpUpperBoundDown(2) <= 1920 || PerpLowerBoundDown(2) <= 2000

                % Downwards Search

                % Making a smaller boxed region to search for the next point
                CoordCompUpperDown = PerpCentroid < PerpUpperBoundDown;
                CoordCompLowerDown = PerpLowerBoundDown < PerpCentroid;
                % Identify corner points for an imaginary scquare around the next point
                UpperPattern_maskDown = (CoordCompUpperDown(:, 1) == 1) & (CoordCompUpperDown(:, 2) == 0);
                LowerPattern_maskDown = (CoordCompLowerDown(:, 1) == 1) & (CoordCompLowerDown(:, 2) == 0);
                % Use logical indexing to keep the point that matches
                PerpNextPtDown = PerpCentroid(UpperPattern_maskDown & LowerPattern_maskDown, :);

                % Checking to see if that Point exists
                if isempty(PerpNextPtDown)
                    PerpNextPtDown = (PerpLowerBoundDown + PerpUpperBoundDown)./2;
                    %plot(PerpNextPtDown(1),PerpNextPtDown(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                    PerpLowerBoundDown = PerpNextPtDown + [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundDown = PerpNextPtDown + [4*DiameterQD,  2*DiameterQD];



                else
                    PerpCounterDown = PerpCounterDown + 1;
                    allPerpNextPtsDown(PerpCounterDown,:) = PerpNextPtDown;
                    PerpLowerBoundDown = PerpNextPtDown + [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundDown = PerpNextPtDown + [4*DiameterQD,  2*DiameterQD];

                end

            end

            while PerpUpperBoundUp(2) >= -100 || PerpLowerBoundUp(2) >= 0

                % Upwards Search

                % Making a smaller boxed region to search for the next point
                CoordCompUpperUp = PerpCentroid < PerpUpperBoundUp;
                CoordCompLowerUp = PerpLowerBoundUp< PerpCentroid;
                % Identify corner points for an imaginary scquare around the next point
                UpperPattern_maskUp = (CoordCompUpperUp(:, 1) == 1) & (CoordCompUpperUp(:, 2) == 0);
                LowerPattern_maskUp = (CoordCompLowerUp(:, 1) == 1) & (CoordCompLowerUp(:, 2) == 0);
                % Use logical indexing to keep the point that matches
                PerpNextPtUp = PerpCentroid(UpperPattern_maskUp & LowerPattern_maskUp, :);

                % Checking to see if that Point exists
                if isempty(PerpNextPtUp)
                    PerpNextPtUp = (PerpLowerBoundUp + PerpUpperBoundUp)./2;
                    PerpLowerBoundUp = PerpNextPtUp - [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundUp = PerpNextPtUp - [4*DiameterQD,  2*DiameterQD];

                else
                    PerpCounterUp = PerpCounterUp + 1;
                    allPerpNextPtsUp(PerpCounterUp,:) = PerpNextPtUp;
                    PerpLowerBoundUp = PerpNextPtUp - [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundUp = PerpNextPtUp - [4*DiameterQD,  2*DiameterQD];
                end
            end


            % Trim the allPerpNextPts array to the actual number of points added
            allPerpPts = vertcat(PerpStartingPt,allPerpNextPtsDown(1:PerpCounterDown, :),allPerpNextPtsUp(1:PerpCounterUp,:));
        end

        function [x_main,y,x_perp,y_perp,red_transpose_x,red_transpose_y,blue_transpose_x,blue_transpose_y,b,perp_b,m,perp_m] = GridLines(obj,allNextPts,allPerpPts,grayImage,scaling,num_lines,sep_red,sep_blue)
            % Description:
                % - creates a line of best fit for each axis based on the QD found through both diagonals 
            % Inputs:
                % AllNextPts -  QD on the diagonal of image for red axis (x-axis)
                % AllPerpPts -  QD on the other diagonal of image for blue axis (y-axis)
                % grayImage - raw image gray scaled 
                % scaling - image compressing factor  
                % num_lines - number of lines parellel to each axis line 
                % sep_red - seperation of red lines (hard coded although can be changed by user)
                % sep_blue - seperation of blue lines (hard coded although can be changed by user)
            % Outputs:
                % All variables return their respective names for each line 

            xvals = allNextPts(:, 1);
            yvals = allNextPts(:, 2);

            % Fitting a line using polyfit
            p = polyfit(xvals, yvals, 1);
            m = p(1);
            b = p(2);

            % Defining the y values of the line
            x = linspace(1, 2560,100);
            x_main = x;
            y = m .* x + b;

            % Trimming outside visible region
            xOOB = x > width(grayImage) / scaling;
            yOOB = (y > height(grayImage) / scaling) | y < 0;
            MainOOB = xOOB | yOOB;
            x_main(MainOOB) = NaN;
            y(MainOOB) = NaN;

            % Generate and plot parallel lines for the red line
            red_all_y_shift = zeros(2*num_lines+1,100);
            red_all_x = repmat(x,2*num_lines+1,1);
            red_line_counter = 0;
            for i = -num_lines:num_lines
                b_shift = b + i * sep_red;
                y_shift = m * x + b_shift;
                red_line_counter = red_line_counter + 1;
                red_all_y_shift(red_line_counter,:) = y_shift;
            end
            red_transpose_x = red_all_x';
            red_transpose_y = red_all_y_shift';

            % Determine out-of-bounds (OOB) points
            xOOB = red_all_x > width(grayImage) / scaling;
            yOOB = (red_all_y_shift > height(grayImage) / scaling) | red_all_y_shift < 0;
            OOB = xOOB | yOOB;

            % Set OOB points to NaN
            red_transpose_x(OOB') = NaN; % Transpose OOB to match the transposed matrices
            red_transpose_y(OOB') = NaN; % Transpose OOB to match the transposed matrices


            % Perpendicular slope
            %
            if length(allPerpPts) < 5
                x_perp = x;
                perp_m = -1 / m;  % Negative reciprocal of the original slope
                length_next_pts = length(allNextPts);
                perp_b = allNextPts(round(length_next_pts./2), 2) - perp_m * allNextPts(round(length_next_pts./2), 1);
                y_perp = perp_m .* x_perp + perp_b;
                fprintf("Using The calculated slope method\n")
            else
                fprintf("Using the QD filtered slope method\n")
                % X and Y values of Second Axes
                Perp_xvals = allPerpPts(:, 1);
                Perp_yvals = allPerpPts(:, 2);

                % Fitting a line using polyfit
                p_perp = polyfit(Perp_xvals, Perp_yvals, 1);
                perp_m = p_perp(1);
                perp_b = p_perp(2);

                % Defining the y values of the line
                x = linspace(1, 2560,100);
                x_perp = x;
                y_perp = perp_m .* x + perp_b;
            end

            % Trimming outside visible region
            xOOB = x_perp > width(grayImage) / scaling;
            yOOB = (y_perp > height(grayImage) / scaling) | y_perp < 0;
            MainOOB = xOOB | yOOB;
            x_perp(MainOOB) = NaN;
            y_perp(MainOOB) = NaN;


            % Generate and plot parallel lines for the blue line
            blue_all_y_shift = zeros(2*num_lines+1,100);
            blue_all_x = repmat(x,2*num_lines+1,1);
            blue_line_counter = 0;
            for i = -num_lines:num_lines
                b_shift = perp_b + i * sep_blue;
                y_shift = perp_m * x + b_shift;
                blue_line_counter = blue_line_counter + 1;
                blue_all_y_shift(blue_line_counter,:) = y_shift;
            end
            blue_transpose_x = blue_all_x';
            blue_transpose_y = blue_all_y_shift';

            % Determine out-of-bounds (OOB) points
            xOOB = blue_all_x > width(grayImage) / scaling;
            yOOB = (blue_all_y_shift > height(grayImage) / scaling) | blue_all_y_shift < 0;
            OOB = xOOB | yOOB;

            % Set OOB points to NaN
            blue_transpose_x(OOB') = NaN; % Transpose OOB to match the transposed matrices
            blue_transpose_y(OOB') = NaN; % Transpose OOB to match the transposed matrices
        end

        function [VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = VirtualQD(obj,num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img)
            % Description:
                % - using the grid pattern previously solved for, does pattern completion to fill in possible missed or not missed QD spots  
            % Inputs:
                % num_lines - number of lines parellel to each axis line 
                % sep_red - seperation of red lines (hard coded although can be changed by user)
                % sep_blue - seperation of blue lines (hard coded although can be changed by user)
                % variables with name variations of b and m are the slopes and intercepts of their respective axes 
            % Outputs:
                % VirtualQDList -  list of all virtual QD
                % FullQDList_sorted - list of both real and virtual QD
                % AllPossibleQDList - list of all possible positions of QD
                % RealQDCentroids - list of real QD
                % centroidx, centroidy - x and y coordinates of real QD respectively 

            % Pre Allocating virtual QD entries and defining a counter
            VirtualQDList = zeros(100, 2);
            AllPossibleQDList = zeros(100, 2);
            DotCounter = 0;
            DotCounter2 = 0;
            for i = -num_lines:num_lines
                for j = -num_lines:num_lines
                    % Calculate new intercepts
                    b1_shift = b + i * sep_red;
                    b2_shift = perp_b + j * sep_blue;

                    % Solve for x and y to find intersection
                    x_intersect = (b2_shift - b1_shift) / (m - perp_m);
                    y_intersect = m * x_intersect + b1_shift;

                    % Plot circle at intersection
                    % Prevention of Real QD and virtual QD overlap
                    DistComp = sqrt((CopyCentroid(:, 1) - x_intersect).^2 + (CopyCentroid(:, 2) - y_intersect).^2);

                    % Prevention of quantum dots outside of visible region being recorded
                    condition = (x_intersect >= 30 && x_intersect <= size(Img, 2)-30) && (y_intersect >= 30 && y_intersect <= size(Img, 1)-30) && all(DistComp > 40);
                    if condition
                        %viscircles([x_intersect, y_intersect], 30, 'Color', [0, 1, 1,0.5], 'LineWidth', 0.5); % Surronding Radius of center point
                        VirtualQD = [x_intersect,y_intersect];
                        DotCounter = DotCounter + 1;
                        VirtualQDList(DotCounter,:) = VirtualQD;
                    end

                    condition2 = (x_intersect >= 30 && x_intersect <= size(Img, 2)-30) && (y_intersect >= 30 && y_intersect <= size(Img, 1)-30);
                    if condition2
                        %viscircles([x_intersect, y_intersect], 30, 'Color', [0, 1, 1,0.5], 'LineWidth', 0.5); % Surronding Radius of center point
                        AllPossibleQD = [x_intersect,y_intersect];
                        DotCounter2 = DotCounter2 + 1;
                        AllPossibleQDList(DotCounter2,:) = AllPossibleQD;
                    end


                end
            end


            % Getting rid of excess entries in the list
            VirtualQDList = VirtualQDList(1:DotCounter,:);
            AllPossibleQDList = AllPossibleQDList(1:DotCounter2,:);



            % Filtering Outlying QD points detected as real
            rows_to_delete_RealQD = false(size(CopyCentroid,1),1);
            for filt = 1:size(CopyCentroid,1)
                DistComp_OutlierQD = sqrt((AllPossibleQDList(:, 1) - CopyCentroid(filt,1)).^2 + (AllPossibleQDList(:, 2) - CopyCentroid(filt,2)).^2);
                minDistOutlierQD = min(DistComp_OutlierQD);
                if minDistOutlierQD >= 50
                    rows_to_delete_RealQD(filt) = true;
                end
            end
            CopyCentroid(rows_to_delete_RealQD,:) = [];

            centroidx = CopyCentroid(:,1);
            centroidy = CopyCentroid(:,2);


            % Filtering Virtual QD points intersecting real QD points
            rows_to_delete_VirtualQD = false(size(VirtualQDList,1),1);
            for filtVirtual = 1:size(VirtualQDList,1)
                DistComp_OutlierVirtualQD = sqrt((CopyCentroid(:, 1) - VirtualQDList(filtVirtual,1)).^2 + (CopyCentroid(:, 2) - VirtualQDList(filtVirtual,2)).^2);
                minDistOutlierVirtualQD = min(DistComp_OutlierVirtualQD);
                if minDistOutlierVirtualQD <= 50 %&& all(DistComp_OutlierQD <=50)
                    rows_to_delete_VirtualQD(filtVirtual)= true;
                end
            end
            VirtualQDList(rows_to_delete_VirtualQD,:) =[];


            % Putting Virtual and Real points into one list
            FullQDList = vertcat(VirtualQDList,CopyCentroid);
            FullQDList_sorted = sortrows(FullQDList, 2);
            RealQDCentroids = CopyCentroid;

        end

        function [DistanceBetweenPoints,NextPt] = RasterScanPatt(obj,radiusQD,StartingQD,FullQDList_sorted,direction,repeatingNum)
             % Description:
                % - path finding for ANC300, user gives movement direction and the path finder finds the next QD in that direction  
            % Inputs:
                % radiusQD - hardcoded approximate radius of each QD
                % startingQD - QD that the path finding starts from
                % FullQDlist_sorted - list of real and virtual QD 
                % direction - direction that the path finder looks for QD. Acceptable inputs are: 'bottomright','bottomleft','topright','topleft'
                % repeatingNum - number of QD found in that direction  
            % Outputs:
                % DistanceBetweenPoints - X and Y distance between each QD within specificed direction
                % Nextpt - centroid coordinate of Next QD within specified direction  

            % For direction can choose bottomright, bottomleft, upperright, upperleft
            DiameterQD = 2*radiusQD;
            DistanceBetweenPoints = zeros(repeatingNum,2);
            counter = 0;

            switch lower(direction)
                case 'bottomright'
                    LowerBound = StartingQD + [2*DiameterQD,  4*DiameterQD];
                    UpperBound = StartingQD + [4*DiameterQD,  2*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 1) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 1) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        LowerBound = NextPt + [2*DiameterQD,  4*DiameterQD];
                        UpperBound = NextPt + [4*DiameterQD,  2*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;
                    end

                case 'bottomleft'
                    LowerBound = StartingQD - [2*DiameterQD, - 4*DiameterQD];
                    UpperBound = StartingQD - [4*DiameterQD, - 2*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        LowerBound = NextPt - [2*DiameterQD, - 4*DiameterQD];
                        UpperBound = NextPt - [4*DiameterQD, - 2*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;

                    end

                case 'topright'
                    LowerBound = StartingQD + [4*DiameterQD, - 2*DiameterQD];
                    UpperBound = StartingQD + [2*DiameterQD, - 4*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);

                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        LowerBound = NextPt + [4*DiameterQD, - 2*DiameterQD];
                        UpperBound = NextPt + [2*DiameterQD, - 4*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;
                    end

                case 'topleft'
                    LowerBound = StartingQD - [4*DiameterQD, 2*DiameterQD];
                    UpperBound = StartingQD - [2*DiameterQD, 4*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 1) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 1) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        LowerBound = NextPt - [4*DiameterQD, 2*DiameterQD];
                        UpperBound = NextPt - [2*DiameterQD, 4*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;
                    end
            end
        end

        function [DistanceBetweenPoints,NextPt,AllPointsOnDiagonal] = Axes_RasterScanPatt(obj,radiusQD,StartingQD,FullQDList_sorted,direction,repeatingNum)
            % Description:
                % - path finding for ANC300, user gives movement direction and the path finder finds the next QD in that direction  (similar to RasterScanPatt with additional features)
            % Inputs:
                % radiusQD - hardcoded approximate radius of each QD
                % startingQD - QD that the path finding starts from
                % FullQDlist_sorted - list of real and virtual QD 
                % direction - direction that the path finder looks for QD. Acceptable inputs are: 'bottomright','bottomleft','topright','topleft'
                % repeatingNum - number of QD found in that direction  
            % Outputs:
                % DistanceBetweenPoints - X and Y distance between each QD within specificed direction
                % Nextpt - centroid coordinate of Next QD within specified direction  
                % AllpointsOnDiagonal - returns every single next pt in order of which they were found in specific direction 
           
           
           
            % For direction can choose bottomright, bottomleft, upperright, upperleft
            DiameterQD = 2*radiusQD;
            DistanceBetweenPoints = zeros(repeatingNum,2);
            counter = 0;
            AllPointsOnDiagonal = [];

            switch lower(direction)
                case 'bottomright'
                    LowerBound = StartingQD + [2*DiameterQD,  4*DiameterQD];
                    UpperBound = StartingQD + [4*DiameterQD,  2*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        %plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        %plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 1) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 1) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        if isempty(NextPt)
                            break
                        end
                        LowerBound = NextPt + [2*DiameterQD,  4*DiameterQD];
                        UpperBound = NextPt + [4*DiameterQD,  2*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal = vertcat(AllPointsOnDiagonal,NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;
                    end

                case 'bottomleft'
                    LowerBound = StartingQD - [2*DiameterQD, - 4*DiameterQD];
                    UpperBound = StartingQD - [4*DiameterQD, - 2*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        %plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        %plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        if isempty(NextPt)
                            break
                        end
                        LowerBound = NextPt - [2*DiameterQD, - 4*DiameterQD];
                        UpperBound = NextPt - [4*DiameterQD, - 2*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal = vertcat(AllPointsOnDiagonal,NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;

                    end

                case 'topright'
                    LowerBound = StartingQD + [4*DiameterQD, - 2*DiameterQD];
                    UpperBound = StartingQD + [2*DiameterQD, - 4*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        %plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        %plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);

                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        if isempty(NextPt)
                            break
                        end
                        LowerBound = NextPt + [4*DiameterQD, - 2*DiameterQD];
                        UpperBound = NextPt + [2*DiameterQD, - 4*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal = vertcat(AllPointsOnDiagonal,NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;
                    end

                case 'topleft'
                    LowerBound = StartingQD - [4*DiameterQD, 2*DiameterQD];
                    UpperBound = StartingQD - [2*DiameterQD, 4*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        counter = counter + 1;
                        %plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        %plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = FullQDList_sorted < UpperBound;
                        CoordCompLower = LowerBound < FullQDList_sorted;
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 1) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 1) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt = FullQDList_sorted(UpperPattern_mask & LowerPattern_mask, :);
                        if isempty(NextPt)
                            break
                        end
                        LowerBound = NextPt - [4*DiameterQD, 2*DiameterQD];
                        UpperBound = NextPt - [2*DiameterQD, 4*DiameterQD];
                        % Finding X and Y distance for ANC300 Controller movement
                        XYdistance = PreviousPt - NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal = vertcat(AllPointsOnDiagonal,NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints(counter,:) = XYdistance;
                    end
            end
            AllPointsOnDiagonal = vertcat(StartingQD,AllPointsOnDiagonal);
        end

        function ANC300Movement(obj,Voltage_X,Frequency_X,Voltage_Y,Frequency_Y,Step_number_X,Step_number_Y,Axis_ID_Direction,ANC300)
            % Description:
                % - sets the voltage and frequency of each axis within piezo and sends movement commands to piezo
            % Inputs:
                % Axis_ID_Direction - stepping direction. Accpetable inputs: "1 up switch", "2 up switch", "1 down switch", "2 down switch"
                % Rest of inputs are self explanatory 
            % Outputs:
                % No output variable as function only writes and sends serial commands to the ANC300   
            
            switch lower(Axis_ID_Direction)
                case "1 up switch"
                    axis = "X";
                    % Setting each axes unused to grounded mode and active axis to stepping mode
                    fprintf(ANC300,"setm 2 gnd");
                    fprintf(ANC300,"setm 3 gnd");
                    fprintf(ANC300,"setm 1 stp");

                    % Setting respective Frequency and Voltage Values
                    Freq_axis_X = "setf 1" + " " + num2str(Frequency_X);
                    fprintf(ANC300,Freq_axis_X);

                    Volt_axis_X = "setv 1" + " " + num2str(Voltage_X);
                    fprintf(ANC300,Volt_axis_X);

                    % Sending Command to start Stepping
                    Step_amount_x = "stepu 1" + " " + num2str(Step_number_X);
                    fprintf("The stepping process is commencing for axis %s\n\n",axis)
                    startPt = tic;
                    fprintf(ANC300,Step_amount_x);
                    fprintf(ANC300,"stepw 1");
                    total_time = toc(startPt);
                    fprintf("Stepping is complete and it took %.5f seconds\n",total_time);

                case "1 down switch"
                    axis = "X";
                    % Setting each axes unused to grounded mode and active axis to stepping mode
                    fprintf(ANC300,"setm 2 gnd");
                    fprintf(ANC300,"setm 3 gnd");
                    fprintf(ANC300,"setm 1 stp");

                    % Setting respective Frequency and Voltage Values
                    Freq_axis_X = "setf 1" + " " + num2str(Frequency_X);
                    fprintf(ANC300,Freq_axis_X);

                    Volt_axis_X = "setv 1" + " " + num2str(Voltage_X);
                    fprintf(ANC300,Volt_axis_X);

                    % Sending Command to start Stepping
                    Step_amount_x = "stepd 1" + " " + num2str(Step_number_X);
                    fprintf("The stepping process is commencing for axis %s\n\n",axis)
                    startPt = tic;
                    fprintf(ANC300,Step_amount_x);
                    fprintf(ANC300,"stepw 1");
                    total_time = toc(startPt);
                    fprintf("Stepping is complete and it took %.5f seconds\n",total_time);

                case "2 up switch"
                    axis = "Y";
                    % Setting each axes unused to grounded mode and active axis to stepping mode
                    fprintf(ANC300,"setm 1 gnd");
                    fprintf(ANC300,"setm 3 gnd");
                    fprintf(ANC300,"setm 2 stp");

                    % Setting respective Frequency and Voltage Values
                    Freq_axis_Y = "setf 2" + " " + num2str(Frequency_Y);
                    fprintf(ANC300,Freq_axis_Y);

                    Volt_axis_Y = "setv 2" + " " + num2str(Voltage_Y);
                    fprintf(ANC300,Volt_axis_Y);

                    % Sending Command to start Stepping
                    Step_amount_Y = "stepu 2" + " " + num2str(Step_number_Y);
                    fprintf("The stepping process is commencing for axis %s\n\n",axis)
                    startPt = tic;
                    fprintf(ANC300,Step_amount_Y);
                    fprintf(ANC300,"stepw 2");
                    total_time = toc(startPt);
                    fprintf("Stepping is complete and it took %.5f seconds\n",total_time);

                case "2 down switch"
                    axis = "Y";
                    % Setting each axes unused to grounded mode and active axis to stepping mode
                    fprintf(ANC300,"setm 1 gnd");
                    fprintf(ANC300,"setm 3 gnd");
                    fprintf(ANC300,"setm 2 stp");

                    % Setting respective Frequency and Voltage Values
                    Freq_axis_Y = "setf 2" + " " + num2str(Frequency_Y);
                    fprintf(ANC300,Freq_axis_Y);

                    Volt_axis_Y = "setv 2" + " " + num2str(Voltage_Y);
                    fprintf(ANC300,Volt_axis_Y);

                    % Sending Command to start Stepping
                    Step_amount_Y = "stepd 2" + " " + num2str(Step_number_Y);
                    fprintf("The stepping process is commencing for axis %s\n\n",axis)
                    startPt = tic;
                    fprintf(ANC300,Step_amount_Y);
                    fprintf(ANC300,"stepw 2");
                    total_time = toc(startPt);
                    fprintf("Stepping is complete and it took %.5f seconds\n",total_time);
                otherwise
                    error('Invalid Axis_ID. Please provide "1" for X-axis or "2" for Y-axis.\n');
            end
        end

        function Dual_ANC300_Movement(obj,Pixel_number_X,Pixel_number_Y,Axis_ID_Direction,ANC300,Frequency,x_factor,y_factor)
             % Description:
                % - similar to ANC300Movement function, however, allows simultaneous movement in both axes without setting voltage and frequency (faster compared to ANC300Movement)
            % Inputs:
                % Axis_ID_Direction - stepping direction. Accpetable inputs: "bottomright", "bottomleft", "topleft", "topright"
                % Rest of inputs are self explanatory 
            % Outputs:
                % No output variable as function only writes and sends serial commands to the ANC300   
            
            
            % Transformation factor from pixel to step
            X_Stepping = abs(Pixel_number_X  * x_factor); % pixels * (steps/pixels)
            Y_Stepping = abs(Pixel_number_Y * y_factor); % pixels * (steps/pixels)
            X_Stepping = round(X_Stepping); % pixels * (steps/pixels)
            Y_Stepping = round(Y_Stepping); % pixels * (steps/pixels)

            % print statements for number of steps (uncomment for debugging
            % purposes) 
            % fprintf("X_step: %d\n",X_Stepping)
            % fprintf("Y_step: %d\n",Y_Stepping)

            if X_Stepping >= 500 | Y_Stepping >= 500 
                error("LARGE STEPPING ERROR, Call Van Damn")
            end


            switch lower(Axis_ID_Direction)
                case 'bottomright'
                    X_Serial_Comd = sprintf("stepu 1 %d",X_Stepping);
                    if X_Stepping ~= 0 
                    fprintf(ANC300,X_Serial_Comd);
                    end
                    StepQueue(obj,X_Stepping,Frequency); 
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"1")
                    Y_Serial_Comd = sprintf("stepu 2 %d",Y_Stepping);
                    if Y_Stepping ~= 0 
                    fprintf(ANC300,Y_Serial_Comd);
                    end                    
                    StepQueue(obj,Y_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"2")


                case 'bottomleft'
                    X_Serial_Comd = sprintf("stepd 1 %d",X_Stepping);
                    if X_Stepping ~= 0 
                    fprintf(ANC300,X_Serial_Comd);
                    end
                    StepQueue(obj,X_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"1")
                    Y_Serial_Comd = sprintf("stepu 2 %d",Y_Stepping);
                    if Y_Stepping ~= 0 
                    fprintf(ANC300,Y_Serial_Comd);
                    end
                    StepQueue(obj,Y_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"2")

                case 'topright'
                    X_Serial_Comd = sprintf("stepu 1 %d",X_Stepping);
                    if X_Stepping ~= 0 
                    fprintf(ANC300,X_Serial_Comd);
                    end
                    StepQueue(obj,X_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"1")
                    Y_Serial_Comd = sprintf("stepd 2 %d",Y_Stepping);
                    if Y_Stepping ~= 0 
                    fprintf(ANC300,Y_Serial_Comd);
                    end
                    StepQueue(obj,Y_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"2")

                case 'topleft'
                    X_Serial_Comd = sprintf("stepd 1 %d",X_Stepping);
                    if X_Stepping ~= 0 
                    fprintf(ANC300,X_Serial_Comd);
                    end
                    StepQueue(obj,X_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"1")
                    Y_Serial_Comd = sprintf("stepd 2 %d",Y_Stepping);
                    if Y_Stepping ~= 0 
                    fprintf(ANC300,Y_Serial_Comd);
                    end
                    StepQueue(obj,Y_Stepping,Frequency);
                    %UpdateText(obj,X_Stepping,Y_Stepping, time_to_pause,"2")
            end
       
        end
       
        function [n, m] = findSubplotDims(obj,numPlots)
            % Find the closest factors of numPlots that are as close to a square as possible
            
            % Initialize the best dimensions
            n = 1;
            m = numPlots;

            % Loop through possible values for n and m
            for i = 1:ceil(sqrt(numPlots))
                if mod(numPlots, i) == 0
                    % If i is a factor, find the corresponding m
                    j = numPlots / i;
                    % Update n and m to minimize the aspect ratio difference
                    if abs(i - j) < abs(n - m)
                        n = i;
                        m = j;
                    end
                end
            end

            % If n*m is less than numPlots, increase n or m
            if n * m < numPlots
                if n < m
                    n = n + 1;
                else
                    m = m + 1;
                end
            end
        end

        % Mainly Used in ANC300_Finalized_Big_Loop
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        function [grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = Finalized_Analyzed_QDBinaryImg_With_Dots(obj,scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,I)
            % Description: takes raw image and identifies Quantum Dots (QD) alongside their centroids via preprocessing and postprocessing 
            % Inputs:
                % Img - the raw image matrix
                % scaling - used for scaling down image for quicker processing 
                % radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area - parameters used for preprocessing and postprocessing ROI 
            % Outputs: 
                % [grayImage,centroid,CopyCentroid,Contrast_adjusted_Img,Img,MaskedImage - matrices that return different stages of the image processing 
                % centroid,copycentroid - return a m x 2 matrix consisting of the x and y positions of the QD centroids 
                % Mask, Mask_scaled - logical matrices used for creating the ROI within the frame 

            Img = im2gray(I);
            grayImage = imresize(Img, scaling);
            sizeImg = size(grayImage);
            sizeImgScaled = size(Img);
         


            % Get coordinates of the circle.
            angles = linspace(0, 2*pi, 10000);
            x = cos(angles) * radii_big_circle + center_big_circle(1);
            y = sin(angles) * radii_big_circle + center_big_circle(2);
            % Get a mask of the circle
            mask = poly2mask(x, y, sizeImg(1), sizeImg(2));

            % Scaled Up mask
            x_scaled = x ./ scaling;
            y_scaled = y ./ scaling;
            mask_scaled = poly2mask(x_scaled,y_scaled,sizeImgScaled(1),sizeImgScaled(2));


            % Making the shadowing across the image more equal
            FlatFieldAdjustedImg = imflatfield(grayImage,sigma_flatfield);

            % Define the radius for the thin lining
            liningRadius = 10;


            % Dilate the mask to create the outer lining
            se = strel('disk', liningRadius);
            dilatedMask = imdilate(mask, se);

            % Get the thin outer lining by subtracting the original mask from the dilated mask
            thinLining = dilatedMask & ~mask;

            % Apply mask to clear outside image area
            maskedImage = FlatFieldAdjustedImg;
            maskedImage(~mask) = 0;

            % Making the constrast between the dots and the background more clear
            Contrast_adjusted_Img = maskedImage;
            Contrast_adjusted_Img(mask == 1) = imadjust(maskedImage(mask == 1));
            Contrast_adjusted_Img(thinLining) = 255;

            % Getting Rid of Salt-Pepper Noise
            imgFiltered = medfilt2(Contrast_adjusted_Img,Salt_pepper_pixel_factor);

            % Binarizing the image very specifically
            BwImg = imbinarize(imgFiltered, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.19);

            % Reapply mask outside of circle
            RectanglularBWImg = BwImg;
            RectanglularBWImg(~mask) = 1;

            % Reverse the colour of foreground and background
            ReversedBWImg = ~RectanglularBWImg;
            %imshow(ReversedBWImg)

            % Removing any points or pixels that are anomolies
            cleanedImg = bwareaopen(ReversedBWImg,Min_circle_area);


            % Taking the Properties of Remaining Image
            [labeledImage, numObjects] = bwlabel(cleanedImg);
            props = regionprops(labeledImage, 'Area', 'Perimeter', 'Eccentricity','EulerNumber');


            % Initialize a mask to keep desired objects
            filteredMask = false(size(cleanedImg));

            % Define thresholds for circularity and eccentricity
            circularityThreshold = 0.5; % Adjust as needed
            eccentricityThreshold = 0.83; % Adjust as needed


            for k = 1 : numObjects
                area = props(k).Area;
                perimeter = props(k).Perimeter;
                eccentricity = props(k).Eccentricity;
                Hollowness = props(k).EulerNumber; 

                % Calculate circularity
                circularity = 4 * pi * area / (perimeter^2);

                % Filter based on circularity and eccentricity
                if circularity > circularityThreshold && eccentricity < eccentricityThreshold && Hollowness == 1
                    filteredMask(labeledImage == k) = true;
                end
            end

            % Further Filtering Using properties
            [labeledImageV2, numObjectsV2] = bwlabel(filteredMask);
            props_centroid = regionprops(filteredMask, 'Centroid','Area'); % Labelling each centroid
            Centroids = vertcat(props_centroid.Centroid); % Putting them into a list for ease of access
            avg_A_detected_spots = sum([props_centroid.Area])/numObjectsV2; 


            % Doing a minimum distance check to filter out points that have shape and
            % intensity properties but don't follow grid pattern
            filterOutlierMask = false(size(filteredMask));
            for j = 1:numObjectsV2
                DistComp_OutlierQD = sqrt((Centroids(:, 1) - Centroids(j,1)).^2 + (Centroids(:, 2) - Centroids(j,2)).^2);
                SortedDistCompQD = sort(DistComp_OutlierQD);
                MinimumDistCompQD = SortedDistCompQD(2);
                if MinimumDistCompQD <= 80
                    filterOutlierMask(labeledImageV2 == j) = false;
                elseif MinimumDistCompQD <= 180
                    filterOutlierMask(labeledImageV2 == j) = true;
                end
                if props_centroid(j).Area < avg_A_detected_spots/2
                   filterOutlierMask(labeledImageV2 == j) = false; 
                end
            end

            % Array of coordinates scaled and scaling back up
            [FinalbwQDImg,~] = bwlabel(filterOutlierMask);
            statsDots = regionprops(FinalbwQDImg, 'Centroid');
            centroid = vertcat(statsDots.Centroid) ./ scaling;
            CopyCentroid = centroid;

            % Masked Image that will be Used for backgrounds
            GrayedImage = im2gray(Img);
            MaskedImage = imadjust(GrayedImage);
            MaskedImage(~mask_scaled) = 0;
        end

        function [allNextPts,allPerpPts] = Finalized_MainAxes(obj,Img,CopyCentroid,radiusQD)
            % Description:
                % - finds a red axis and a blue axis that act as the corresponding x and y axis for the grid pattern of the QD
            % Inputs:
                % Img - raw image matrix
                % Copycentroid - centroid coordinates of the Quantum Dots
                % radiusQD - hardcoded radius of QD for plotting 
            % Outputs:
                % AllNextPts - returns QD on the diagonal of image for red axis (x-axis)
                % AllPerpPts - returns QD on the other diagonal of image for blue axis (y-axis)

            % Calculate distances from each point to the top right corner
            distancesX = sqrt((CopyCentroid(:, 1) - size(Img,2)).^2);
            distancesY = sqrt((CopyCentroid(:, 2) - 0).^2);
            distances = distancesX + distancesY;

            % Creating a table with the indexs of each point alongside the distance from the top corner
            Indices = 1:length(distances);
            Table_with_dist_ind = table(Indices(:),distances(:),'VariableNames',["Indices", "Distance_to_Corner"]);
            Sorted_Table_with_dist_ind = sortrows(Table_with_dist_ind,"Distance_to_Corner");

            % Estimate the maximum number of points (initial size of coords as an estimate)
            maxPoints = size(CopyCentroid, 1);

            % Parameter for Accuracy Wanted
            num_of_QD_diagonal = 5;

            i = 0;
            while i <= length(distances)
                i = i + 1;
                % Get the closest point coordinates
                closestPoint = CopyCentroid(Sorted_Table_with_dist_ind{i,1}, :);
                StartingPt = closestPoint;

                % Finding next closest point
                closestPointNew = closestPoint + [radiusQD, -radiusQD];

                % Preallocate array to store all next points
                allNextPts = zeros(maxPoints+1, 2);


                % For Box pattern creating a lower boundary to prevent skewing of Axes line
                DiameterQD = 2*radiusQD;
                LowerBound = closestPointNew - [2.5*DiameterQD, - 4.5*DiameterQD];
                UpperBound = closestPointNew - [4.5*DiameterQD, - 2.5*DiameterQD];

                % Redefining each set variable
                centroid = CopyCentroid;
                numNextPts = 0;


                while ~isempty(centroid)
                    % Remove points equal to the closest point
                    targetRow = all(centroid == closestPoint, 2);
                    centroid(targetRow, :) = [];

                    % Check if coords is empty after removal
                    if isempty(centroid)
                        break;
                    end

                    % Making a smaller boxed region below the previous point to search for the next point
                    CoordCompUpper = centroid < UpperBound;
                    CoordCompLower = LowerBound < centroid;
                    % Identify corner points for an imaginary scquare around the next point
                    UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                    LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);

                    % Use logical indexing to keep the point that matches
                    NextPt = centroid(UpperPattern_mask & LowerPattern_mask, :);
                    centroid(UpperPattern_mask & LowerPattern_mask, :) = [];

                    % Add the closest point to the list of all next points
                    if ~isempty(NextPt)
                        numNextPts = numNextPts + 1;
                        CopynumNextPts = numNextPts;
                        allNextPts(numNextPts, :) = NextPt;
                    end

                    % Use logical indexing to remove all the points before each lowerbound
                    % to prevent the while loop from going on forever
                    RightOfLowerBound = LowerBound < centroid;
                    PointsRemovedRightSide = (RightOfLowerBound(:,1) == 1);
                    centroid(PointsRemovedRightSide,:) = [];

                    % Use logical indexing to remove all the points above each upperbound
                    % to prevent the while loop from going on forever
                    AboveUpperBound = UpperBound > centroid;
                    PointsRemovedAbove = (AboveUpperBound(:,2) == 1);
                    centroid(PointsRemovedAbove,:) = [];

                    % Check if coords is empty after filtering
                    if isempty(centroid)
                        break;
                    end


                    if isempty(NextPt)
                        LowerBound = closestPointNew - [6*DiameterQD, - 8*DiameterQD];
                        UpperBound = closestPointNew - [8*DiameterQD, - 5*DiameterQD];
                        closestPointNew = (UpperBound + LowerBound)./2;
                    else
                        % Update the closestPoint, closestPointNew, Lowerbound for the next iteration
                        closestPoint = NextPt;
                        closestPointNew = closestPoint + [radiusQD, -radiusQD];
                        LowerBound = closestPointNew - [2.5*DiameterQD, - 4.5*DiameterQD];
                        UpperBound = closestPointNew - [4.5*DiameterQD, - 2.5*DiameterQD];
                    end


                end
                if numNextPts + 1 >= num_of_QD_diagonal
                    break
                elseif (i == length(distances)) && (numNextPts < num_of_QD_diagonal)
                    num_of_QD_diagonal = num_of_QD_diagonal - 1;
                    i = 0;
                elseif num_of_QD_diagonal == 0
                    error("No possible QD were found on the diagonal")
                end

            end

            % Trim the allNextPts array to the actual number of points added
            allNextPts = vertcat(StartingPt,allNextPts(1:CopynumNextPts, :));

            % Find QD across the Perpendicular diagonal
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            PerpCentroid = CopyCentroid;

            PerpStartingPt = allNextPts(round(length(allNextPts)/2),:);
            PerpLowerBoundDown =  PerpStartingPt + [2*DiameterQD,  4*DiameterQD];
            PerpUpperBoundDown = PerpStartingPt + [4*DiameterQD,  2*DiameterQD];

            PerpLowerBoundUp = PerpStartingPt - [4*DiameterQD, 2*DiameterQD];
            PerpUpperBoundUp = PerpStartingPt - [2*DiameterQD, 4*DiameterQD];


            % Preallocate array to store all next points
            allPerpNextPtsDown = zeros(maxPoints, 2);
            allPerpNextPtsUp = zeros(maxPoints, 2);

            % Counter to trim list
            PerpCounterDown = 0;
            PerpCounterUp = 0;

            while PerpUpperBoundDown(2) <= 1920 || PerpLowerBoundDown(2) <= 2000

                % Making a smaller boxed region to search for the next point
                CoordCompUpperDown = PerpCentroid < PerpUpperBoundDown;
                CoordCompLowerDown = PerpLowerBoundDown < PerpCentroid;
                % Identify corner points for an imaginary scquare around the next point
                UpperPattern_maskDown = (CoordCompUpperDown(:, 1) == 1) & (CoordCompUpperDown(:, 2) == 0);
                LowerPattern_maskDown = (CoordCompLowerDown(:, 1) == 1) & (CoordCompLowerDown(:, 2) == 0);
                % Use logical indexing to keep the point that matches
                PerpNextPtDown = PerpCentroid(UpperPattern_maskDown & LowerPattern_maskDown, :);

                % Checking to see if that Point exists
                if isempty(PerpNextPtDown)
                    PerpNextPtDown = (PerpLowerBoundDown + PerpUpperBoundDown)./2;
                    PerpLowerBoundDown = PerpNextPtDown + [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundDown = PerpNextPtDown + [4*DiameterQD,  2*DiameterQD];


                else
                    PerpCounterDown = PerpCounterDown + 1;
                    allPerpNextPtsDown(PerpCounterDown,:) = PerpNextPtDown;
                    PerpLowerBoundDown = PerpNextPtDown + [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundDown = PerpNextPtDown + [4*DiameterQD,  2*DiameterQD];
                end

            end

            while PerpUpperBoundUp(2) >= -100 || PerpLowerBoundUp(2) >= 0

                % Upwards Search
                %p5 = plot(PerpLowerBoundUp(1),PerpLowerBoundUp(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[1, 0, 1],"LineStyle","none");
                %p6 = plot(PerpUpperBoundUp(1),PerpUpperBoundUp(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 1, 0],"LineStyle","none");
                % Making a smaller boxed region to search for the next point
                CoordCompUpperUp = PerpCentroid < PerpUpperBoundUp;
                CoordCompLowerUp = PerpLowerBoundUp< PerpCentroid;
                % Identify corner points for an imaginary scquare around the next point
                UpperPattern_maskUp = (CoordCompUpperUp(:, 1) == 1) & (CoordCompUpperUp(:, 2) == 0);
                LowerPattern_maskUp = (CoordCompLowerUp(:, 1) == 1) & (CoordCompLowerUp(:, 2) == 0);
                % Use logical indexing to keep the point that matches
                PerpNextPtUp = PerpCentroid(UpperPattern_maskUp & LowerPattern_maskUp, :);

                % Checking to see if that Point exists
                if isempty(PerpNextPtUp)
                    PerpNextPtUp = (PerpLowerBoundUp + PerpUpperBoundUp)./2;
                    %plot(PerpNextPtUp(1),PerpNextPtUp(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                    PerpLowerBoundUp = PerpNextPtUp - [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundUp = PerpNextPtUp - [4*DiameterQD,  2*DiameterQD];

                else
                    PerpCounterUp = PerpCounterUp + 1;
                    allPerpNextPtsUp(PerpCounterUp,:) = PerpNextPtUp;
                    PerpLowerBoundUp = PerpNextPtUp - [2*DiameterQD,  4*DiameterQD];
                    PerpUpperBoundUp = PerpNextPtUp - [4*DiameterQD,  2*DiameterQD];
                    %plot(PerpNextPtUp(1),PerpNextPtUp(2),"Marker","o","MarkerSize",5,"MarkerFaceColor","auto")
                end
            end


            % Trim the allPerpNextPts array to the actual number of points added
            allPerpPts = vertcat(PerpStartingPt,allPerpNextPtsDown(1:PerpCounterDown, :),allPerpNextPtsUp(1:PerpCounterUp,:));
        end

        function [x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = Finalized_GridLines(obj,allNextPts,allPerpPts)
             % Description:
                % - creates a line of best fit for each axis based on the QD found through both diagonals 
            % Inputs:
                % AllNextPts -  QD on the diagonal of image for red axis (x-axis)
                % AllPerpPts -  QD on the other diagonal of image for blue axis (y-axis)
            % Outputs:
                % All variables return their respective name properties for each line
            
            
            % Primary Axis
            xvals = allNextPts(:, 1);
            yvals = allNextPts(:, 2);

            % Fitting a line using polyfit
            p = polyfit(xvals, yvals, 1);
            m = p(1);
            b = p(2);

            % Defining the y values of the line
            x = linspace(1, 2560,100);
            x_main = x;
            y = m .* x + b;


            % Perpendicular Axis to primary axis
            % If enough points don't exist manually calculate a line
            if length(allPerpPts) < 5
                x_perp = x;
                perp_m = -1 / m;  % Negative reciprocal of the original slope
                length_next_pts = length(allNextPts);
                perp_b = allNextPts(round(length_next_pts./2), 2) - perp_m * allNextPts(round(length_next_pts./2), 1);
                y_perp = perp_m .* x_perp + perp_b;
                %fprintf("Using The calculated slope method\n")

                % If enough points exist use the ones filtered diagonally
            else
                x_perp = [];
                y_perp = []; 
                %fprintf("Using the QD filtered slope method\n")
                % X and Y values of Second Axes
                Perp_xvals = allPerpPts(:, 1);
                Perp_yvals = allPerpPts(:, 2);

                % Fitting a line using polyfit
                p_perp = polyfit(Perp_xvals, Perp_yvals, 1);
                perp_m = p_perp(1);
                perp_b = p_perp(2);
            end
        end

        function [VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = Finalized_VirtualQD(obj,num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled)
            % Description:
                % - using the grid pattern previously solved for, does pattern completion to fill in missed or not missed QD spots not occupied by real QD 
            % Inputs:
                % num_lines - number of lines parellel to each axis line 
                % sep_red - seperation of red lines (hard coded although can be changed by user)
                % sep_blue - seperation of blue lines (hard coded although can be changed by user)
                % variables with name variations of b and m are the slopes and intercepts of their respective axes 
            % Outputs:
                % VirtualQDList -  list of all virtual QD
                % FullQDList_sorted - list of both real and virtual QD
                % AllPossibleQDList - list of all possible positions of QD
                % RealQDCentroids - list of real QD
                % centroidx, centroidy - x and y coordinates of real QD respectively 

            % Pre Allocating virtual QD entries and defining a counter
            VirtualQDList = zeros(100, 2);
            AllPossibleQDList = zeros(100, 2);
            DotCounter = 0;
            DotCounter2 = 0;
            for i = -num_lines:num_lines
                for j = -num_lines:num_lines
                    % Calculate new intercepts
                    b1_shift = b + i * sep_red;
                    b2_shift = perp_b + j * sep_blue;

                    % Solve for x and y to find intersection
                    x_intersect = (b2_shift - b1_shift) / (m - perp_m);
                    y_intersect = m * x_intersect + b1_shift;

                    % Plot circle at intersection
                    % Prevention of Real QD and virtual QD overlap
                    DistComp = sqrt((CopyCentroid(:, 1) - x_intersect).^2 + (CopyCentroid(:, 2) - y_intersect).^2);

                    % Prevention of quantum dots outside of visible region being recorded
                    condition = (x_intersect >= 1 && x_intersect <= size(Img, 2)-1) && (y_intersect >= 1 && y_intersect <= size(Img, 1)-1) && all(DistComp > 40);
                    if condition
                        if mask_scaled(ceil(y_intersect),ceil(x_intersect))
                        %viscircles([x_intersect, y_intersect], 30, 'Color', [0, 1, 1,0.5], 'LineWidth', 0.5); % Surronding Radius of center point
                        VirtualQD = [x_intersect,y_intersect];
                        DotCounter = DotCounter + 1;
                        VirtualQDList(DotCounter,:) = VirtualQD;
                        end
                    end

                    condition2 = (x_intersect >= 1 && x_intersect <= size(Img, 2)-1) && (y_intersect >= 1 && y_intersect <= size(Img, 1)-1);
                    if condition2
                        if mask_scaled(ceil(y_intersect),ceil(x_intersect))
                        %viscircles([x_intersect, y_intersect], 30, 'Color', [0, 1, 1,0.5], 'LineWidth', 0.5); % Surronding Radius of center point
                        AllPossibleQD = [x_intersect,y_intersect];
                        DotCounter2 = DotCounter2 + 1;
                        AllPossibleQDList(DotCounter2,:) = AllPossibleQD;
                        end
                    end


                end
            end


            % Getting rid of excess entries in the list
            VirtualQDList = VirtualQDList(1:DotCounter,:);
            AllPossibleQDList = AllPossibleQDList(1:DotCounter2,:);



            % Filtering Outlying QD points detected as real
            rows_to_delete_RealQD = false(size(CopyCentroid,1),1);
            for filt = 1:size(CopyCentroid,1)
                DistComp_OutlierQD = sqrt((AllPossibleQDList(:, 1) - CopyCentroid(filt,1)).^2 + (AllPossibleQDList(:, 2) - CopyCentroid(filt,2)).^2);
                minDistOutlierQD = min(DistComp_OutlierQD);
                if minDistOutlierQD >= 50
                    rows_to_delete_RealQD(filt) = true;
                end
            end
            CopyCentroid(rows_to_delete_RealQD,:) = [];

            centroidx = CopyCentroid(:,1);
            centroidy = CopyCentroid(:,2);


            % Filtering Virtual QD points intersecting real QD points
            rows_to_delete_VirtualQD = false(size(VirtualQDList,1),1);
            for filtVirtual = 1:size(VirtualQDList,1)
                DistComp_OutlierVirtualQD = sqrt((CopyCentroid(:, 1) - VirtualQDList(filtVirtual,1)).^2 + (CopyCentroid(:, 2) - VirtualQDList(filtVirtual,2)).^2);
                minDistOutlierVirtualQD = min(DistComp_OutlierVirtualQD);
                if minDistOutlierVirtualQD <= 50 %&& all(DistComp_OutlierQD <=50)
                    rows_to_delete_VirtualQD(filtVirtual)= true;
                end
            end
            VirtualQDList(rows_to_delete_VirtualQD,:) =[];


            % Putting Virtual and Real points into one list
            RealQDCentroids = CopyCentroid;
            FullQDList = vertcat(VirtualQDList,RealQDCentroids);
            FullQDList_sorted = sortrows(FullQDList, 2);

        end

        % this is the old non optimized raster scan function not in use but don't delete in case of debugging! 
        function [DistanceBetweenPoints_Original,DistanceBetweenPoints_Rotated,NextPt_Coords,Rotated_NextPt,AllPointsOnDiagonal_Original,AllPointsOnDiagonal_Rotated] = Finalized_RasterScanPatt(obj,radiusQD,StartingQD,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted,direction,repeatingNum)
            % Description:
                % - path finding for ANC300, user gives movement direction and the path finder finds the next QD in that direction  
            % Inputs:
                % radiusQD - hardcoded approximate radius of each QD
                % startingQD - QD that the path finding starts from
                % Rotated_Table_FullQDList_sorted - rotated centroid list of real and virtual QD that match movement axes of ANC300 
                % Table_FullQDList_sorted - centoid list of all real and virtual QD 
                % direction - direction that the path finder looks for QD. Acceptable inputs are: 'bottomright','bottomleft','topright','topleft'
                % repeatingNum - number of QD found in that direction  
            % Outputs:
                % DistanceBetweenPoints - X and Y distance between each QD within specificed direction
                % Nextpt - centroid coordinate of Next QD within specified direction  
                % AllpointsOnDiagonal - returns every single next pt in order of which they were found in specific direction 
                % All named outputs have their rotated counterparts as outputs as well 
           
            % Defining variables and creating lists 
            DiameterQD = 2*radiusQD;
            DistanceBetweenPoints_Original = zeros(repeatingNum,2);
            counter = 0;
            AllPointsOnDiagonal_Original = [];
            AllPointsOnDiagonal_Rotated = [];
            [~,~,ID_StartingQD] = intersect(StartingQD,Table_FullQDList_sorted.("QD Coordinates"),"rows");
            rotated_StartingQD = Rotated_Table_FullQDList_sorted{ID_StartingQD,"QD Coordinates"};


            switch lower(direction)
                case 'bottomright'
                    LowerBound = StartingQD + [2*DiameterQD,  6*DiameterQD];
                    UpperBound = StartingQD + [6*DiameterQD,  2*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        [~,~,ID_PrevPt] = intersect(StartingQD,Table_FullQDList_sorted.("QD Coordinates"),"rows");
                        counter = counter + 1;
                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = Table_FullQDList_sorted.("QD Coordinates") < UpperBound;
                        CoordCompLower = LowerBound < Table_FullQDList_sorted.("QD Coordinates");
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 1) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 1) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt_Mask = UpperPattern_mask & LowerPattern_mask;
                        NextPt_Coords = Table_FullQDList_sorted{NextPt_Mask,"QD Coordinates"};
                        NextPt_ID = Table_FullQDList_sorted{NextPt_Mask,"ID"};
                        LowerBound = NextPt_Coords + [2*DiameterQD,  4*DiameterQD];
                        UpperBound = NextPt_Coords + [4*DiameterQD,  2*DiameterQD];
                        % Finding X and Y distance of original frame
                        XYdistance_original = PreviousPt - NextPt_Coords;
                        % Finding X and Y distance of rotated frame
                        Rotated_PreviousPt = Rotated_Table_FullQDList_sorted{ID_PrevPt,"QD Coordinates"};
                        Rotated_NextPt = Rotated_Table_FullQDList_sorted{NextPt_ID,"QD Coordinates"};
                        XYdistance_rotated = Rotated_PreviousPt - Rotated_NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt_Coords;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal_Original = vertcat(AllPointsOnDiagonal_Original,NextPt_Coords);
                        AllPointsOnDiagonal_Rotated = vertcat(AllPointsOnDiagonal_Rotated,Rotated_NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints_Original(counter,:) = XYdistance_original;
                        DistanceBetweenPoints_Rotated(counter,:) = XYdistance_rotated;
                    end

                case 'bottomleft'
                    LowerBound = StartingQD - [2*DiameterQD, - 6*DiameterQD];
                    UpperBound = StartingQD - [6*DiameterQD, - 2*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum
                        [~,~,ID_PrevPt] = intersect(StartingQD,Table_FullQDList_sorted.("QD Coordinates"),"rows");
                        counter = counter + 1;
                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);
                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = Table_FullQDList_sorted.("QD Coordinates") < UpperBound;
                        CoordCompLower = LowerBound < Table_FullQDList_sorted.("QD Coordinates");
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt_Mask = UpperPattern_mask & LowerPattern_mask;
                        NextPt_Coords = Table_FullQDList_sorted{NextPt_Mask,"QD Coordinates"};
                        NextPt_ID = Table_FullQDList_sorted{NextPt_Mask,"ID"};

                        % Logical boundery finding
                        LowerBound = NextPt_Coords - [2*DiameterQD, - 6*DiameterQD];
                        UpperBound = NextPt_Coords - [6*DiameterQD, - 2*DiameterQD];
                        % Finding X and Y distance of original frame
                        XYdistance_original = PreviousPt - NextPt_Coords;
                        % Finding X and Y distance of rotated frame
                        Rotated_PreviousPt = Rotated_Table_FullQDList_sorted{ID_PrevPt,"QD Coordinates"};
                        Rotated_NextPt = Rotated_Table_FullQDList_sorted{NextPt_ID,"QD Coordinates"};
                        XYdistance_rotated = Rotated_PreviousPt - Rotated_NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        %PreviousPt = NextPt_Coords;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal_Original = vertcat(AllPointsOnDiagonal_Original,NextPt_Coords);
                        AllPointsOnDiagonal_Rotated = vertcat(AllPointsOnDiagonal_Rotated,Rotated_NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints_Original(counter,:) = XYdistance_original;
                        DistanceBetweenPoints_Rotated(counter,:) = XYdistance_rotated;
                    end

                case 'topright'
                    LowerBound = StartingQD + [6*DiameterQD, - 2*DiameterQD];
                    UpperBound = StartingQD + [2*DiameterQD, - 6*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum

                        [~,~,ID_PrevPt] = intersect(StartingQD,Table_FullQDList_sorted.("QD Coordinates"),"rows");
                        counter = counter + 1;

                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);

                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = Table_FullQDList_sorted.("QD Coordinates") < UpperBound;
                        CoordCompLower = LowerBound < Table_FullQDList_sorted.("QD Coordinates");
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 0) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 0) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt_Mask = UpperPattern_mask & LowerPattern_mask;
                        NextPt_Coords = Table_FullQDList_sorted{NextPt_Mask,"QD Coordinates"};
                        NextPt_ID = Table_FullQDList_sorted{NextPt_Mask,"ID"};

                        LowerBound = NextPt_Coords + [4*DiameterQD, - 2*DiameterQD];
                        UpperBound = NextPt_Coords + [2*DiameterQD, - 4*DiameterQD];

                        % Finding X and Y distance of original frame
                        XYdistance_original = PreviousPt - NextPt_Coords;
                        % Finding X and Y distance of rotated frame
                        Rotated_PreviousPt = Rotated_Table_FullQDList_sorted{ID_PrevPt,"QD Coordinates"};
                        Rotated_NextPt = Rotated_Table_FullQDList_sorted{NextPt_ID,"QD Coordinates"};
                        XYdistance_rotated = Rotated_PreviousPt - Rotated_NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt_Coords;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal_Original = vertcat(AllPointsOnDiagonal_Original,NextPt_Coords);
                        AllPointsOnDiagonal_Rotated = vertcat(AllPointsOnDiagonal_Rotated,Rotated_NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints_Original(counter,:) = XYdistance_original;
                        DistanceBetweenPoints_Rotated(counter,:) = XYdistance_rotated;
                    end

                case 'topleft'
                    LowerBound = StartingQD - [6*DiameterQD, 2*DiameterQD];
                    UpperBound = StartingQD - [2*DiameterQD, 6*DiameterQD];
                    PreviousPt = StartingQD;
                    while counter < repeatingNum

                        [~,~,ID_PrevPt] = intersect(StartingQD,Table_FullQDList_sorted.("QD Coordinates"),"rows");
                        counter = counter + 1;

                        plot(LowerBound(1),LowerBound(2),"Marker","o","LineWidth",2,"MarkerSize",5,"MarkerEdgeColor",[0.2010 0.7450 0.9330]);
                        plot(UpperBound(1),UpperBound(2),"Marker","o","LineWidth",2,"Markersize",5,"MarkerEdgeColor",[1, 0.75, 0.2]);

                        % Making a smaller boxed region to search for the next point
                        CoordCompUpper = Table_FullQDList_sorted.("QD Coordinates") < UpperBound;
                        CoordCompLower = LowerBound < Table_FullQDList_sorted.("QD Coordinates");
                        % Identify corner points for an imaginary scquare around the next point
                        UpperPattern_mask = (CoordCompUpper(:, 1) == 1) & (CoordCompUpper(:, 2) == 0);
                        LowerPattern_mask = (CoordCompLower(:, 1) == 1) & (CoordCompLower(:, 2) == 0);
                        % Use logical indexing to keep the point that matches
                        NextPt_Mask = UpperPattern_mask & LowerPattern_mask;
                        NextPt_Coords = Table_FullQDList_sorted{NextPt_Mask,"QD Coordinates"};
                        NextPt_ID = Table_FullQDList_sorted{NextPt_Mask,"ID"};

                        LowerBound = NextPt_Coords - [4*DiameterQD, 2*DiameterQD];
                        UpperBound = NextPt_Coords - [2*DiameterQD, 4*DiameterQD];
                        % Finding X and Y distance of original frame
                        XYdistance_original = PreviousPt - NextPt_Coords;
                        % Finding X and Y distance of rotated frame
                        Rotated_PreviousPt = Rotated_Table_FullQDList_sorted{ID_PrevPt,"QD Coordinates"};
                        Rotated_NextPt = Rotated_Table_FullQDList_sorted{NextPt_ID,"QD Coordinates"};
                        XYdistance_rotated = Rotated_PreviousPt - Rotated_NextPt;
                        % Sketching Raster Scan Path
                        %plot(vertcat(PreviousPt(:,1),NextPt(:,1)), vertcat(PreviousPt(:,2),NextPt(:,2)), 'LineStyle', '-', 'Color', [0.4, 0, 0.6,0.4], 'LineWidth',2.5)
                        PreviousPt = NextPt_Coords;
                        % Add NextPt to a list of values
                        AllPointsOnDiagonal_Original = vertcat(AllPointsOnDiagonal_Original,NextPt_Coords);
                        AllPointsOnDiagonal_Rotated = vertcat(AllPointsOnDiagonal_Rotated,Rotated_NextPt);
                        % Add to list for ease of use later
                        DistanceBetweenPoints_Original(counter,:) = XYdistance_original;
                        DistanceBetweenPoints_Rotated(counter,:) = XYdistance_rotated;
                    end
            end
            AllPointsOnDiagonal_Original = vertcat(StartingQD,AllPointsOnDiagonal_Original);
            AllPointsOnDiagonal_Rotated = vertcat(rotated_StartingQD,AllPointsOnDiagonal_Rotated);



        end

        function [NextPoint,Rotated_NextPoint,XY_Difference,Rotated_XY_Difference] = OptimizedRasterScan(obj,startPoint, Centroids_Table, Rotated_Centroids_Table,direction)
            % Description:
                % - path finding for ANC300, user gives movement direction and the path finder finds the next QD in that direction  
            % Inputs:
                % Startpoint - QD that the path finding starts from
                % Rotated_Centroids_Table - rotated centroid list of real and virtual QD that match movement axes of ANC300 
                % Centroids_Table - centoid list of all real and virtual QD 
                % direction - direction that the path finder looks for QD. Acceptable inputs are: 'bottomright','bottomleft','topright','topleft'
            % Outputs:
                % Rotated_NextPoint - rotated centroid coordinate of Next QD within specified direction
                % Nextpt - centroid coordinate of Next QD within specified direction  
                % XY_Difference - X and Y pixel difference between starting point and next point 
                % Rotated_XY_Difference - rotated X and Y pixel difference between starting point and next point 

            % defining rotated points table and other need variables
            [~,~,ID_StartingQD] = intersect(startPoint,Centroids_Table.("QD Coordinates"),"rows");
            rotated_StartingPoint = Rotated_Centroids_Table{ID_StartingQD,"QD Coordinates"};
            
            Centroids = Centroids_Table.("QD Coordinates");
            Rotated_Centroids = Rotated_Centroids_Table.("QD Coordinates"); 


            % Extract the numeric data from the table
            x = Centroids(:, 1);
            y = Centroids(:, 2);

            % Exclude the starting point from the list of points
            isDifferentFromStart = (x ~= startPoint(1)) | (y ~= startPoint(2));

            % Apply the condition for the selected direction
            switch direction
                case 'bottomleft'
                    % Points with smallest x and largest y
                    validIndices = (x < startPoint(1)) & (y > startPoint(2)) & isDifferentFromStart;
                    
                case 'bottomright'
                    % Points with largest x and largest y
                    validIndices = (x > startPoint(1)) & (y > startPoint(2)) & isDifferentFromStart;
                    
                case 'topleft'
                    % Points with smallest x and smallest y
                    validIndices = (x < startPoint(1)) & (y < startPoint(2)) & isDifferentFromStart;
                    
                case 'topright'
                    % Points with largest x and smallest y
                    validIndices = (x > startPoint(1)) & (y < startPoint(2)) & isDifferentFromStart;
                    
                otherwise
                    error('Invalid direction input');
            end
            
            % Calculate distances for the valid points
            distances = (x(validIndices) - startPoint(1)).^2 + (y(validIndices) - startPoint(2)).^2;
            validPoints = Centroids(validIndices, :);
            
            if isempty(validPoints)
                error('No valid points found in the specified direction');
            end
            
            % Find the closest point among the valid points (this is the next point now) 
            [~, minIndex] = min(distances);
            NextPoint = validPoints(minIndex, :);
            [~,~,ID_NextPt] = intersect(NextPoint,Centroids,"rows");
            Rotated_NextPoint = Rotated_Centroids(ID_NextPt,:); 

            % calculating the x and y distance between the original and rotated
            % pathway 
            XY_Difference = startPoint - NextPoint; 
            Rotated_XY_Difference = rotated_StartingPoint - Rotated_NextPoint; 
    
        end
        
        function [StartingQD,StartingQD_rotated,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = Precision_Locking(obj,ANC300,PhotoType,QD_counter,pyueye_initialization_return,accuracy_margin)
            % Takes all previous functions that do the image analysis, pattern completion, and movement and put it to a more organized formating running a constnat while loop looking for the dot we need 

            % Determining if the photo being taken is the starting QD determiner or stepping QD determiner 
            if PhotoType == "StartingPhoto"
                ImageName = "StartingPhoto.jpg"; % Photo File name
                %pyrun_file_text = sprintf("simple_pyueye_snap.py '%s'",ImageName); % Photo Snapping File name used for legacy snapping 
                title_text = sprintf('QD Starting Point: %s\n',ImageName); % Figure title name 
            elseif PhotoType == "SteppingPhoto"
                ImageName = sprintf("Nanowire_Photo_[%d %d].jpg",QD_counter); % Photo File name
                %pyrun_file_text = sprintf("simple_pyueye_snap.py '%s'",ImageName);  % Photo Snapping File name used for legacy snapping 
                title_text = sprintf('QD Stepping Point: %s\n',ImageName); % Figure title name 
            else
                error("Invalid Photo Type entered: Try again\n")
            end

            % Filtering Settings
            scaling = 0.5;
            skyBlue = [0.53, 0.81, 0.92];
            radii_big_circle = 500;
            center_big_circle = [650,500];
            sigma_flatfield = 35;
            Salt_pepper_pixel_factor = [20 20];
            Min_circle_area = 200;
            radiusQD = 25; %Pre Aug 9th 2024,30;


            % Hard coded location of LED Spot
            LEDSpotCentroid = LED_Cooordinate_Identifer(obj,"","Read","LAB"); 
            LEDSpotCentroid = round(LEDSpotCentroid); 

            % Define the number of parallel lines and their separation
            num_lines = 9; % Num lines on each side of the original line
            sep_red = 345; % Separation distance red lines
            sep_blue = 370; % Seperation distance blue lines

            % other paramters of use
            Frequency = 20; 
            angle = 45.1; % Pre Aug 16th 43.696474776653320; % use the angle found by the algorithm in the XYANC300Axes script 

            % X and y factors
            Read_XY_factor  = XY_Factor_Identifier(obj,"","Read","LAB"); 
            x_factor = Read_XY_factor.X_factor; 
            y_factor = Read_XY_factor.Y_factor; 
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            Eucli_ShortestDistance = accuracy_margin*2; % Variable set for purposes of initializing while loop 
            Iterations = 0; 


            
            while Eucli_ShortestDistance > accuracy_margin 
            Iterations = Iterations + 1; 

            % Snapping Photo 
            py.pyueye_func.snap_image(pyueye_initialization_return, ImageName)
            %Fetch today's folder string 
            date = py.qd_data_folder_creation.date_string(); 
            date = string(date);
            date_path = sprintf("C:\\Users\\Quantum Dot\\Desktop\\Bera Yavuz - ANC300 Movement and Images\\QD_Data\\%s_Test\\Position_uEYE",date);
            addpath(date_path)
            I = ImageName; 
        

            % Identifying real QD
            [grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = Finalized_Analyzed_QDBinaryImg_With_Dots(obj,scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,I);
            
            % Finds QD across the diagonal of image and collects data to form two main axes
            [allNextPts,allPerpPts] = Finalized_MainAxes(obj,Img,CopyCentroid,radiusQD);
            % Creates grid lines to find intersections
            [x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = Finalized_GridLines(obj,allNextPts,allPerpPts); 
            % Creates and plots Virtual and Real QD
            [VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = Finalized_VirtualQD(obj,num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled); 
            % Rotates Img and dots displayed 
            [VirtualQDList_rotated,RealQD_rotated,LEDSpotCentroid_rotated,Table_FullQDList_sorted,Rotated_Table_FullQDList_sorted,rotated_image] = Rotated_Img_n_Pts(obj,angle,MaskedImage,FullQDList_sorted,VirtualQDList,RealQDCentroids,LEDSpotCentroid);
            %rotated_image = imrotate(MaskedImage, angle,"bilinear","crop");
            % Finding the Starting Point based on Image 
            [StartingQD_rotated,ShortestDistance,ID] = FindClosestPt(obj,Rotated_Table_FullQDList_sorted.("QD Coordinates"),LEDSpotCentroid_rotated); 
            Eucli_ShortestDistance = EucliDistance(obj,LEDSpotCentroid_rotated,StartingQD_rotated); 
            StartingQD = Table_FullQDList_sorted.("QD Coordinates")(ID,:); 
            %fprintf("QD Coord Difference: X = %.3f | Y = %.3f\n",ShortestDistance)
            %fprintf("Closest QD distance: %.2f\n",Eucli_ShortestDistance)

            % figure("Name",title_text,"Color",skyBlue)
            % imshow(rotated_image)
            % hold on; 
            % title_text_iter = sprintf(title_text + "iteration: %d",Iterations); 
            % title(title_text_iter)
            % plot(VirtualQDList_rotated(:,1), VirtualQDList_rotated(:,2), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r',"LineStyle","none"); % Rotated Virtual QD plotted
            % plot(RealQD_rotated(:,1), RealQD_rotated(:,2), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none"); % Rotated Real QD plotted
            % plot(LEDSpotCentroid_rotated(1),LEDSpotCentroid_rotated(2),"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o","LineStyle","none"); % Rotated LED plotted 
            % plot(StartingQD_rotated(1),StartingQD_rotated(2), 'bo', 'MarkerSize', 10,'MarkerFaceColor','b','LineStyle',"none");
            % legend("Virtual QD", "Real QD", "LED spot", "Starting QD")
            % 

            % Finding direction need to travel to in order to get to QD 
            [direction] = ClosestPtDirection(obj,ShortestDistance); 
            %fprintf("%s\n",direction)
            % Moving to the startingQD
            Dual_ANC300_Movement(obj,ShortestDistance(1),ShortestDistance(2),direction,ANC300,Frequency,x_factor,y_factor)
            if Iterations > 50
                %fprintf("%---------------------------------------\n%d iterations were done\n",Iterations)
                break
            end
            end 
            
            
            %fprintf("%---------------------------------------\n%d iterations were done\n",Iterations)
            % fprintf("total iterations took %.2f seconds\n---------------------------------------\n",elaseped_Time)
            % 

 
            
        end
 
        function[StartingQD,StartingQD_rotated,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = Precision_Locking_Matlab(obj,ANC300,QD_counter,vid_UI,src_UI,accuracy_margin)
            % Takes all previous functions that do the image analysis, pattern completion, and movement and put it to a more organized formating running a constnat while loop looking for the dot we need 


            % Filtering Settings
            scaling = 0.5;
            skyBlue = [0.53, 0.81, 0.92];
            radii_big_circle = 500;
            center_big_circle = [650,500];
            sigma_flatfield = 35;
            Salt_pepper_pixel_factor = [20 20];
            Min_circle_area = 200;
            radiusQD = 25; %Pre Aug 9th 2024,30;


            % Hard coded location of LED Spot
            LEDSpotCentroid = LED_Cooordinate_Identifer(obj,"","Read","LAB"); 
            LEDSpotCentroid = round(LEDSpotCentroid); 

            % Define the number of parallel lines and their separation
            num_lines = 9; % Num lines on each side of the original line
            sep_red = 345; % Separation distance red lines
            sep_blue = 370; % Seperation distance blue lines

            % other paramters of use
            Frequency = 20; 
            angle = 45.1; % Pre Aug 16th 43.696474776653320; % use the angle found by the algorithm in the XYANC300Axes script 

            % X and y factors
            Read_XY_factor  = XY_Factor_Identifier(obj,"","Read","LAB"); 
            x_factor = Read_XY_factor.X_factor; 
            y_factor = Read_XY_factor.Y_factor; 
            %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            Eucli_ShortestDistance = accuracy_margin*2; % Variable set for purposes of initializing while loop 
            Iterations = 0; 


            
            while Eucli_ShortestDistance > accuracy_margin
            Iterations = Iterations + 1; 

            % Snapping Photo 
            
            [UI_Position_Img] = UI_Snap_Img(obj,vid_UI,src_UI,"No",QD_counter); 
        

            % Identifying real QD
            [grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = Finalized_Analyzed_QDBinaryImg_With_Dots(obj,scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,UI_Position_Img);
            
            % Finds QD across the diagonal of image and collects data to form two main axes
            [allNextPts,allPerpPts] = Finalized_MainAxes(obj,Img,CopyCentroid,radiusQD);
            % Creates grid lines to find intersections
            [x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = Finalized_GridLines(obj,allNextPts,allPerpPts); 
            % Creates and plots Virtual and Real QD
            [VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = Finalized_VirtualQD(obj,num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled); 
            % Rotates Img and dots displayed 
            [VirtualQDList_rotated,RealQD_rotated,LEDSpotCentroid_rotated,Table_FullQDList_sorted,Rotated_Table_FullQDList_sorted,rotated_image] = Rotated_Img_n_Pts(obj,angle,MaskedImage,FullQDList_sorted,VirtualQDList,RealQDCentroids,LEDSpotCentroid);
            %rotated_image = imrotate(MaskedImage, angle,"bilinear","crop");
            % Finding the Starting Point based on Image 
            [StartingQD_rotated,ShortestDistance,ID] = FindClosestPt(obj,Rotated_Table_FullQDList_sorted.("QD Coordinates"),LEDSpotCentroid_rotated); 
            Eucli_ShortestDistance = EucliDistance(obj,LEDSpotCentroid_rotated,StartingQD_rotated); 
            StartingQD = Table_FullQDList_sorted.("QD Coordinates")(ID,:); 
            %fprintf("QD Coord Difference: X = %.3f | Y = %.3f\n",ShortestDistance)
            %fprintf("Closest QD distance: %.2f\n",Eucli_ShortestDistance)

            % figure("Name",title_text,"Color",skyBlue)
            % imshow(rotated_image)
            % hold on; 
            % title_text_iter = sprintf(title_text + "iteration: %d",Iterations); 
            % title(title_text_iter)
            % plot(VirtualQDList_rotated(:,1), VirtualQDList_rotated(:,2), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r',"LineStyle","none"); % Rotated Virtual QD plotted
            % plot(RealQD_rotated(:,1), RealQD_rotated(:,2), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none"); % Rotated Real QD plotted
            % plot(LEDSpotCentroid_rotated(1),LEDSpotCentroid_rotated(2),"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o","LineStyle","none"); % Rotated LED plotted 
            % plot(StartingQD_rotated(1),StartingQD_rotated(2), 'bo', 'MarkerSize', 10,'MarkerFaceColor','b','LineStyle',"none");
            % legend("Virtual QD", "Real QD", "LED spot", "Starting QD")
            % 

            % Finding direction need to travel to in order to get to QD 
            [direction] = ClosestPtDirection(obj,ShortestDistance); 
            %fprintf("%s\n",direction)
            % Moving to the startingQD
            Dual_ANC300_Movement(obj,ShortestDistance(1),ShortestDistance(2),direction,ANC300,Frequency,x_factor,y_factor)
            pause(0.2)
            if Iterations > 50
                %fprintf("%---------------------------------------\n%d iterations were done\n",Iterations)
                break
            end
            end 
            
            
            %fprintf("%---------------------------------------\n%d iterations were done\n",Iterations)
            % fprintf("total iterations took %.2f seconds\n---------------------------------------\n",elaseped_Time)
            % 

 
            
        end
 

        % Side Line functions 
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        function [VirtualQDList_rotated,RealQD_rotated,LEDSpotCentroid_rotated,Table_FullQDList_sorted,Rotated_Table_FullQDList_sorted,rotated_image] = Rotated_Img_n_Pts(obj,angle,MaskedImage,FullQDList_sorted,VirtualQDList,RealQDCentroids,LEDSpotCentroid)
            % Description:
                % - rotates everything by a given angle   
            % Inputs:
                % angle - angle of rotation to apply for all things 
                % FullQDList_sorted - centoid list of all real and virtual QD 
                % VirtualQDList - centoid list of all virtal QD 
                % RealQDCentroids - centoid list of all real QD 
                % LEDSpotCentroid - centroid of excitation laser 
                % MaskedImage - raw image with ROI highlighted and constrat adjusted   
            % Outputs:
                % All input centroids rotated 
                % Table_FullQDList_sorted - Organized table with centroid list for real and virtual QD
                % Rotated_Table_FullQDList_sorted - Organized rotated table with centroid list for real and virtual QD
                % rotated_image - raw image rotated 


            % Creating a table with all the points and assigning them an ID Val
            num_IDS = length(FullQDList_sorted);
            Table_FullQDList_sorted = table((1:num_IDS)',FullQDList_sorted,'VariableNames',{'ID','QD Coordinates'});

            % Rotating the Image
            rotated_image = imrotate(MaskedImage, angle,"bilinear","crop");
            center_pt = [size(MaskedImage, 2)/2, size(MaskedImage, 1)/2]; % Image center
            %fprintf("Center Point %f\n",center_pt)

            % Rotate points
            theta = -angle * pi / 180; % Convert angle to radians

            % Virtual points rotated
            VirtualQDList_X_rotated = center_pt(1) + (VirtualQDList(:,1) - center_pt(1)) * cos(theta) - (VirtualQDList(:,2) - center_pt(2)) * sin(theta);
            VirtualQDList_Y_rotated = center_pt(2) + (VirtualQDList(:,1) - center_pt(1)) * sin(theta) + (VirtualQDList(:,2) - center_pt(2)) * cos(theta);
            VirtualQDList_rotated = horzcat(VirtualQDList_X_rotated,VirtualQDList_Y_rotated); 

            % Real points rotated
            RealQD_X_rotated = center_pt(1) + (RealQDCentroids(:,1) - center_pt(1)) * cos(theta) - (RealQDCentroids(:,2) - center_pt(2)) * sin(theta);
            RealQD_Y_rotated = center_pt(2) + (RealQDCentroids(:,1) - center_pt(1)) * sin(theta) + (RealQDCentroids(:,2) - center_pt(2)) * cos(theta);
            RealQD_rotated = horzcat(RealQD_X_rotated,RealQD_Y_rotated);

            % Excitation laser rotated
            LEDSpotCentroidX_rotated = center_pt(1) + (LEDSpotCentroid(1) - center_pt(1)) * cos(theta) - (LEDSpotCentroid(2) - center_pt(2)) * sin(theta);
            LEDSpotCentroidY_rotated = center_pt(2) + (LEDSpotCentroid(1) - center_pt(1)) * sin(theta) + (LEDSpotCentroid(2) - center_pt(2)) * cos(theta);
            LEDSpotCentroid_rotated = horzcat(LEDSpotCentroidX_rotated,LEDSpotCentroidY_rotated);

            % Rotating table that has all points combined 
            Table_FullQDList_sorted_X = Table_FullQDList_sorted.("QD Coordinates")(:,1);
            Table_FullQDList_sorted_Y = Table_FullQDList_sorted.("QD Coordinates")(:,2);
            Rotated_Table_FullQDList_sorted_X = center_pt(1) + (Table_FullQDList_sorted_X - center_pt(1)) * cos(theta) - (Table_FullQDList_sorted_Y - center_pt(2)) * sin(theta);
            Rotated_Table_FullQDList_sorted_Y = center_pt(2) + (Table_FullQDList_sorted_X - center_pt(1)) * sin(theta) + (Table_FullQDList_sorted_Y - center_pt(2)) * cos(theta);
            FullQDList_rotated_coords = horzcat(Rotated_Table_FullQDList_sorted_X,Rotated_Table_FullQDList_sorted_Y);
            Rotated_Table_FullQDList_sorted = table((1:num_IDS)',FullQDList_rotated_coords,'VariableNames',{'ID','QD Coordinates'});
        end
   
        function [ClosestPt,ShortestDistance,ID] = FindClosestPt(obj, AllPossibleQDList,LEDSpotCentroid)
             % Description:
                % - finds closest QD to excitation laser    
            % Inputs:
                % AllPossibleQDList - centroid list of real and virtual QD
                % LEDSpotCentroid - coordinates of excitation laser
            % Outputs:
                % ClosestPt - closest QD to excitation laser 
            
            DistComp = sqrt((AllPossibleQDList(:, 1) - LEDSpotCentroid(1)).^2 + (AllPossibleQDList(:, 2) - LEDSpotCentroid(2)).^2);
            [~,ID] = min(DistComp);
            ShortestDistance = LEDSpotCentroid - AllPossibleQDList(ID,:); 
            ClosestPt = AllPossibleQDList(ID,:); 
        end

        function time_to_pause = StepQueue(obj,Step_num,frequency)
         % Description:
                % - acts as a queue for stepping to prevent the next line of code from being read before stepping is complete 
            % Inputs:
                % Step_num - number of serial steps sent to ANC300
                % frequency - frquency at which the ANC300 is moving 
            % Outputs:
                % No output varialbe as sole purpose of function is to act as a small delay 

            Error_Margin_Factor = 0.3; % in terms of seconds so change as appropriately 
            time_to_pause = (Step_num/frequency) + Error_Margin_Factor; % calculating the time needed for pausing between lines 
            pause(time_to_pause); 
        end

        function [direction] = ClosestPtDirection(obj,ShortestDistance)
            % determines the direction of the closest point to the LED spot for starting QD purposes
            % Inputs:
                % Distance = distance between closest point and LED Spot 
            % Outputs:
                % returns the direction of the closest point. Possible outputs are; "bottomleft","bottomright","topleft","topright" 

                if (ShortestDistance(1) < 0 &&  ShortestDistance(2) < 0)
                    direction = "bottomright"; 
                elseif (ShortestDistance(1) > 0 &&  ShortestDistance(2) < 0)
                    direction = "bottomleft"; 
                elseif (ShortestDistance(1) > 0 &&  ShortestDistance(2) > 0)
                    direction = "topleft"; 
                elseif (ShortestDistance(1) < 0 && ShortestDistance(2) > 0)
                    direction = "topright";
                end
        end
                    
        function UpdateText(obj,X_step_num,Y_step_num, time_to_pause,axis)

            switch axis 
                case "1"
                fprintf("---------------------------------------\nStepping Complete\n%d step(s) in X direction\nExpected time of Completion: %.3f seconds\n---------------------------------------\n",X_step_num,time_to_pause)
                case "2"
                fprintf("---------------------------------------\nStepping Complete\n%d step(s) in Y direction\nExpected time of Completion: %.3f seconds\n---------------------------------------\n",Y_step_num,time_to_pause)    
            end
        end

        function distance = EucliDistance(obj,Point_1,Point_2)
            % Finds the euclidean distance (shortest distance) between two
            % points. 

            distance = sqrt( (Point_1(1)-Point_2(1))^2 + (Point_1(2)-Point_2(2))^2 );
        end

        function readQD_position = QD_tracking_N_Identification(obj,QD_position,InputSource,Read_Write,Device,direction_QD)
            % Description:
                % - % keeps track of the current QD position by keeping track within a text file  
            % Inputs:
                % QD_position - current position of QD in the form of [row,column]
                % InputSource - specifies if the input was automatically put or if the user manually put QD position (valid inputs: "Auto" or "Manual") 
                % Read_Write - specifies if user is trying to get the current position or if current position is being updated (valid inputs: "Read" or "Write")
                % Device - specifies what device user is using in order to access specific QD history file (valid inputs: "HOME_PC", "LAB")
                % direction_QD - the direction at which the raster scan is travelling 
            % Outputs:
                % readQD_position -  current QD position 
            
            % Define the filename
            %filenameHome = 'QD_History.txt';
            if Device == "LAB"
            filenameLAB = "C:\Users\Quantum Dot\Desktop\Bera_Yavuz_GitHub\AttoCube-Project-Stuff\QD_History.txt";
            elseif Device == "HOME_PC"
            filenameLAB = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\QD_History.txt"; 
            elseif Device == "MAC"
                filenameLAB = "/Users/bera_yavuz/Desktop/GitHub ANC300 Project/QD_History.txt"; 
            else
                fprintf("invalid Device name try again")
            end
            % Open the file for appending
            
            % Convert QD position to a string to put in textfile
            QD_position = mat2str(QD_position);

            switch Read_Write 

                case "Write" 
                fileID = fopen(filenameLAB, 'a+');
                % Check if the file opened successfully
                if fileID == -1
                    error('Failed to open the file.');
                end

                % Get the current timestamp
                currentDateTime = datetime('now');
                currentDateTimeText = string(currentDateTime);
                if InputSource == "Auto";
                    InputSourceText = sprintf("QD position was updated --Automatically--");
                elseif InputSource == "Manual"
                    InputSourceText = sprintf("QD position was updated --Manually--");
                else
                    error("Invalid Input: Please Input Auto or Manual")
                end 
                % Write the timestamp and information to the file
                fprintf(fileID, '%s\n================================================\n%s\nRaster Scan Direction: %s\nCurrent position of QD:\n%s\n\n\n', currentDateTimeText, InputSourceText,direction_QD,QD_position);

                
                case "Read"
                fileID = fopen(filenameLAB, 'rt');
                if fileID == -1
                    error('Failed to open the file.');
                end
                
                % Read the entire file into a single string
                fileContent = fread(fileID, '*char')';
                
                % Split the content into lines
                lines = strsplit(fileContent, newline);
                
                % Initialize variables to store the last three non-empty lines
                lastLine = '';
                secondLastLine = '';
                thirdLastLine = '';
                
                nonEmptyCount = 0;  % Counter for non-empty lines
                
                % Loop through the lines in reverse order
                for i = length(lines):-1:1
                    if ~isempty(strtrim(lines{i}))
                        nonEmptyCount = nonEmptyCount + 1;
                        if nonEmptyCount == 1
                            lastLine = lines{i};
                        elseif nonEmptyCount == 2
                            secondLastLine = lines{i};
                        elseif nonEmptyCount == 3
                            thirdLastLine = lines{i};
                            % Split the string into words
                            words = strsplit(thirdLastLine);
                            % Extract the last word
                            Raster_DirectionQD = words{end};
                            break;  % Stop after finding the third-to-last non-empty line
                        end
                    end
                end
                
                % Convert the last and third-to-last lines to numbers
                readQD_position = struct('lastLine', str2num(lastLine), 'thirdLastLine', Raster_DirectionQD); 
                   
            end
            if Read_Write == "Write"
                readQD_position = ''; 
            end
            fclose(fileID);
        end
        
        function read_LED_Coordinate  = LED_Cooordinate_Identifer(obj,LED_coordinate,Read_Write,Device)
            % Description:
                % finds and keeps track of the current LED coordinate 
            % Inputs:
                % LED_coordinate - Specifies x and y coordinate of LED 
                % Read_Write - specifies if user is trying to get the current position or if current position is being updated (valid inputs: "Read" or "Write")
                % Device - specifies what device user is using in order to access specific QD history file (valid inputs: "HOME_PC", "LAB")
            % Outputs:
                % read_LED_Coordinate -  current LED coordinate 
            
            % Define the filename 
            if Device == "LAB"
            filename = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images\LED_Coordinate_History.txt";
            elseif Device == "HOME_PC"
            filename = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\LED_Coordinate_History.txt"; 
            else
                fprintf("invalid Device name try again")
            end
            % Open the file for appending
            
            % Convert QD position to a string to put in textfile
            LED_coordinate = mat2str(LED_coordinate);

            switch Read_Write 

                case "Write" 
                fileID = fopen(filename, 'a+');
                % Check if the file opened successfully
                if fileID == -1
                    error('Failed to open the file.');
                end

                % Get the current timestamp
                currentDateTime = datetime('now');
                currentDateTimeText = string(currentDateTime);
               
                % Write the timestamp and information to the file
                fprintf(fileID, '%s\n================================================\nCurrent LED coordinate:\n%s\n\n\n', currentDateTimeText,LED_coordinate);

                
                case "Read"
                fileID = fopen(filename, 'rt');
                if fileID == -1
                    error('Failed to open the file.');
                end
                
                % Read the entire file into a single string
                fileContent = fread(fileID, '*char')';
                
                % Split the content into lines
                lines = strsplit(fileContent, newline);
                
                % Initialize variables to store the last three non-empty lines
                lastLine = '';

                
                nonEmptyCount = 0;  % Counter for non-empty lines
                
                % Loop through the lines in reverse order
                for i = length(lines):-1:1
                    if ~isempty(strtrim(lines{i}))
                        nonEmptyCount = nonEmptyCount + 1;
                        if nonEmptyCount == 1
                            lastLine = lines{i};
                        end
                    end
                end
                
                % Convert the last and third-to-last lines to numbers
                read_LED_Coordinate = str2num(lastLine); 
                   
            end
            if Read_Write == "Write"
                read_LED_Coordinate = ''; 
            end
            fclose(fileID);
        end
        
        function Read_XY_factor  = XY_Factor_Identifier(obj,XY_factor,Read_Write,Device)
            % Description:
                % finds and keeps track of the current XY factor 
            % Inputs:
                % XY_factor - Specifies x and y factor (steps/pixels unit)  
                % Read_Write - specifies if user is trying to get the current position or if current position is being updated (valid inputs: "Read" or "Write")
                % Device - specifies what device user is using in order to access specific QD history file (valid inputs: "HOME_PC", "LAB")
            % Outputs:
                % Read_XY_factor -  current XY factor coordinate 
            
            % Define the filename
            if Device == "LAB" 
            filename = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\XY_Factor_Images\XY_Factor_History.txt";
            elseif Device == "HOME_PC"
            filename = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\XY_Factor_History.txt"; 
            else
                fprintf("invalid Device name try again")
            end
    

            switch Read_Write 

                case "Write" 
                fileID = fopen(filename, 'a+');
                % Check if the file opened successfully
                if fileID == -1
                    error('Failed to open the file.');
                end

                % Get the current timestamp
                currentDateTime = datetime('now');
                currentDateTimeText = string(currentDateTime);
               
                % Write the timestamp and information to the file
                fprintf(fileID, '%s\n================================================\nCurrent X Factor: %s\nCurrent Y Factor: %s\n\n\n', currentDateTimeText,XY_factor);

                
                case "Read"
                fileID = fopen(filename, 'rt');
                if fileID == -1
                    error('Failed to open the file.');
                end
                
                % Read the entire file into a single string
                fileContent = fread(fileID, '*char')';
                
                % Split the content into lines
                lines = strsplit(fileContent, newline);
                
                % Initialize variables to store the last three non-empty lines
                lastLine = '';
                secondLastLine = '';
               
                
                nonEmptyCount = 0;  % Counter for non-empty lines
                
                % Loop through the lines in reverse order
                for i = length(lines):-1:1
                    if ~isempty(strtrim(lines{i}))
                        nonEmptyCount = nonEmptyCount + 1;
                        if nonEmptyCount == 1
                            lastLine = lines{i};
                        elseif nonEmptyCount == 2
                            secondLastLine = lines{i};
                              % Stop after finding the second-to-last non-empty line
                        end
                    end
                    % Split the string into words
                    words = strsplit(lastLine);
                    % Extract the last word
                    Y_Factor = words{end};
                    % Split the string into words
                    words = strsplit(secondLastLine);
                    % Extract the last word
                    X_Factor = words{end};
                end
                
                % Convert the last and third-to-last lines to numbers
                Read_XY_factor = struct("X_factor",str2num(X_Factor),"Y_factor",str2num(Y_Factor)); 
                   
            end
            if Read_Write == "Write"
                Read_XY_factor = ''; 
            end
            fclose(fileID);
        end

        function Read_Last_Log  = Last_Log_Identification(obj,Read_Write,Device)
            % Description:
                % finds and keeps track of the current and previous QD history file currently in use  
            % Inputs:
                % Read_Write - specifies if user is trying to get the current position or if current position is being updated (valid inputs: "Read" or "Write")
                % Device - specifies what device user is using in order to access specific QD history file (valid inputs: "HOME_PC", "LAB")
            % Outputs:
                % Read_Last_Log -  current QD history text file iteration 
            
            % Define the filename
            if Device == "LAB"
            filename =  "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_History_Logs\QD_History_Log.txt";
            elseif Device == "MAC"
            filename = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\QD_History_log.txt";
            elseif Device == "HOME_PC"
            filename = "C:\Users\yavub\OneDrive\Desktop\QD_History_logs\QD_History_log.txt"; 
            else
                fprintf("invalid Device name try again")
            end
 
            % Get the current timestamp
            currentDateTime = datetime('now');
            currentDateTimeText = string(currentDateTime); % returns accurate year_month_day and time 


            t = datetime("now");
            [y,m,d] = ymd(t);
            date_today = sprintf("_%d_%d_%d",y,m,d);
            HistoryTextFileName = sprintf("QD_History%s.txt",date_today);
    
 
            switch Read_Write
 
                case "Write"
                    fileID = fopen(filename, 'rt');
                    % Check if the file opened successfully
                    if fileID == -1
                        error('Failed to open the file.');
                    end
                    
                    % Read the entire file into a single string
                    fileContent = fread(fileID, '*char')';
   
                    % Split the content into lines
                    lines = strsplit(fileContent, newline);
                    
                    % Initialize variables to store the last three non-empty lines
                    lastLine = '';
                
                    
                    nonEmptyCount = 0;  % Counter for non-empty lines
                    
                    % Loop through the lines in reverse order
                    for i = length(lines):-1:1
                        if ~isempty(strtrim(lines{i}))
                            nonEmptyCount = nonEmptyCount + 1;
                            if nonEmptyCount == 1
                                lastLine = lines{i};
                                break 
                            end
                        end                    
                    end
                    if lastLine == HistoryTextFileName
                        fprintf("Log already exists")
                        fclose(fileID);
                        return
                    end
                    fclose(fileID);
                    fileID = fopen(filename, 'a+');
                    % Check if the file opened successfully
                    if fileID == -1
                        error('Failed to open the file.');
                    end
                    % Write the timestamp and information to the file
                    fprintf(fileID, '%s\n================================================\n%s\n\n', currentDateTimeText,HistoryTextFileName);
            
 
                case "Read"
                    fileID = fopen(filename, 'rt');
                    % Check if the file opened successfully
                    if fileID == -1
                        error('Failed to open the file.');
                    end
                    
                    % Read the entire file into a single string
                    fileContent = fread(fileID, '*char')';
 
                    if isempty(fileContent)
                        error("Nothing is inside file to be read")
                    end
                    
                    % Split the content into lines
                    lines = strsplit(fileContent, newline);
                    
                    % Initialize variables to store the last three non-empty lines
                    lastline = '';
                
                    
                    nonEmptyCount = 0;  % Counter for non-empty lines
                    
                    % Loop through the lines in reverse order
                    for i = length(lines):-1:1
                        if ~isempty(strtrim(lines{i}))
                            nonEmptyCount = nonEmptyCount + 1;
                            if nonEmptyCount == 1
                                lastline = lines{i};
                                break
                            end
                        end                    
                    end
                    Read_Last_Log = lastline;
            end
            fclose(fileID);
        end

        function readQD_position = QD_tracking_N_IdentificationVer2(obj,QD_position,InputSource,Action,Device,direction_QD)
            % Description:
                % - % keeps track of the current QD position by keeping track within a text file  
            % Inputs:
                % QD_position - current position of QD in the form of [row,column]
                % InputSource - specifies if the input was automatically put or if the user manually put QD position (valid inputs: "Auto" or "Manual") 
                % Action - specifies if user is trying to get the current position or if current position is being updated and if new file is being made(valid inputs: "Read","Write","Initialize")
                % Device - specifies what device user is using in order to access specific QD history file (valid inputs: "HOME_PC", "LAB")
                % direction_QD - the direction at which the raster scan is travelling 
            % Outputs:
                % readQD_position -  current QD position 
            
            % Define the filename
            t = datetime("now");
            [y,m,d] = ymd(t);
            date_today = sprintf("_%d_%d_%d",y,m,d);
            HistoryTextFileName = sprintf("QD_History%s.txt",date_today);
            currentTimeStr = datestr(t,'HH:MM:SS'); % only returns hour:minute:second 


            if Device == "LAB"
            filenameLAB = sprintf("C:\\Users\\Quantum Dot\\Desktop\\Bera Yavuz - ANC300 Movement and Images\\QD_History_Logs\\QD_History_Text_Files\\%s",HistoryTextFileName); 
            elseif Device == "HOME_PC"
            filenameLAB = sprintf("C:\\Users\\yavub\\OneDrive\\Desktop\\QD_History_logs\\All_QD_History_Files\\%s",HistoryTextFileName); 
            elseif Device == "MAC"
            filenameLAB = "/Users/bera_yavuz/Desktop/GitHub ANC300 Project/QD_History.txt"; 
            else
                fprintf("invalid Device name try again")
            end
            % Open the file for appending
            
            % Convert QD position to a string to put in textfile
            QD_position = mat2str(QD_position);

            % Reading off of a previously existing text file for first dot 


            switch Action 

                case "Initialize"
                    if exist(filenameLAB, 'file')
                        fprintf('File already exists.');
                        return
                    else
                        fileID_QD = fopen(filenameLAB, 'w+'); % Create the file
                        if fileID_QD == -1
                            error('Failed to create the file.');
                        end
                        fclose(fileID_QD);
                    end
                    fileID_QD = fopen(filenameLAB,'rt'); % opening file in order to read it line by line 
      
                    % Read the entire file into a single string
                    fileContent = fread(fileID_QD, '*char')';
                    
                    if isempty(fileContent)
                        fclose(fileID_QD); 
                        read_QD_History_log = Last_Log_Identification(obj,"Read","LAB"); 
                        filename_QD_History = sprintf("C:\\Users\\Quantum Dot\\Desktop\\Bera Yavuz - ANC300 Movement and Images\\QD_History_Logs\\QD_History_Text_Files\\%s",read_QD_History_log); 
                        %filename_QD_History = sprintf("C:\\Users\\yavub\\OneDrive\\Desktop\\QD_History_logs\\All_QD_History_Files\\%s",read_QD_History_log); 
                        fileID_Log = fopen(filename_QD_History, 'rt'); 
                        if fileID_Log == -1
                            error('Failed to open the file.');
                        end
                        
                        % Read the entire file into a single string
                        fileContent = fread(fileID_Log, '*char')';
                        
                        % Split the content into lines
                        lines = strsplit(fileContent, newline);
                        
                        % Initialize variables to store the last three non-empty lines
                        lastline = '';
                        nonEmptyCount = 0;  % Counter for non-empty lines
                    
                        % Loop through the lines in reverse order
                        for i = length(lines):-1:1
                            if ~isempty(strtrim(lines{i}))
                                nonEmptyCount = nonEmptyCount + 1;
                                if nonEmptyCount == 1
                                    lastline = lines{i};
                                elseif nonEmptyCount == 3
                                    thirdlastline = lines{i}; 
                                    % split the string into words
                                    words = strsplit(thirdlastline);
                                    % extract the last word (direction)
                                    Raster_DirectionQD = words{end}; 
                                    break
                                end
                            end
                        end
                        QD_position_initial = lastline; 
                        fclose(fileID_Log); 
    
                        % opening back the newly created file to add that line
                        fileID = fopen(filenameLAB,'a+');
                        if fileID == -1
                            error('Failed to open the file.');
                        end
                         % Write the timestamp and information to the file
                        fprintf(fileID, '%s\n================================================\nTaken from last QD text file\nRaster Scan Direction: %s\nCurrent position of QD:\n%s\n\n\n', currentTimeStr,Raster_DirectionQD,QD_position_initial);
                        fclose(fileID); 
                        Last_Log_Identification(obj,"Write","LAB"); 
                    end

                case "Write" 
                
                fileID = fopen(filenameLAB, 'a+');
                % Check if the file opened successfully
                if fileID == -1
                    error('Failed to open the file.');
                end
    
                
                if InputSource == "Auto"
                    InputSourceText = sprintf("QD position was updated --Automatically--");
                elseif InputSource == "Manual"
                    InputSourceText = sprintf("QD position was updated --Manually--");
                else
                    error("Invalid Input: Please Input Auto or Manual")
                end 
                % Write the timestamp and information to the file
                fprintf(fileID, '%s\n================================================\n%s\nRaster Scan Direction: %s\nCurrent position of QD:\n%s\n\n\n', currentTimeStr, InputSourceText,direction_QD,QD_position);
                fclose(fileID); 
                
                case "Read"
                fileID = fopen(filenameLAB, 'rt');
                if fileID == -1
                    error('Failed to open the file.');
                end
                
                % Read the entire file into a single string
                fileContent = fread(fileID, '*char')';

                if isempty(fileContent)
                    fprintf("there is nothing in the file, please fix issue")
                    return 
                end
                
                % Split the content into lines
                lines = strsplit(fileContent, newline);
                
                % Initialize variables to store the last three non-empty lines
                lastLine = '';
                secondLastLine = '';
                thirdLastLine = '';
                
                nonEmptyCount = 0;  % Counter for non-empty lines
                
                % Loop through the lines in reverse order
                for i = length(lines):-1:1
                    if ~isempty(strtrim(lines{i}))
                        nonEmptyCount = nonEmptyCount + 1;
                        if nonEmptyCount == 1
                            lastLine = lines{i};
                        elseif nonEmptyCount == 2
                            secondLastLine = lines{i};
                        elseif nonEmptyCount == 3
                            thirdLastLine = lines{i};
                            % Split the string into words
                            words = strsplit(thirdLastLine);
                            % Extract the last word
                            Raster_DirectionQD = words{end};
                            break;  % Stop after finding the third-to-last non-empty line
                        end
                    end
                end
                
                % Convert the last and third-to-last lines to numbers
                readQD_position = struct('lastLine', str2num(lastLine), 'thirdLastLine', Raster_DirectionQD); 
                fclose(fileID); 
            end
            if Action == "Write"
                readQD_position = ''; 
            end
        end

        function read_fast_movement_setting = fast_movement_settings(obj,setting_updated,setting_value,Device, InputSource)
        % Description:
                % - % keeps track of the current QD position by keeping track within a text file  
            % Inputs:
                % setting_updated - which setting the user wants to chance
                % setting_value - what the selected setting is being set to
                % InputSource - specifies if user is trying to get the current position or if current position is being updated and if new file is being made(valid inputs: "Read","Write","Initialize")
                % Device - specifies what device user is using in order to access specific QD history file (valid inputs: "HOME_PC", "LAB")
            % Outputs:
                % read_fast_movement_setting -  reads setting for each respective thing 
        
            % Define the filename
            t = datetime("now");
            [y,m,d] = ymd(t);
            date_today = sprintf("_%d_%d_%d",y,m,d);
            HistoryTextFileName = sprintf("QD_History%s.txt",date_today);
            currentTimeStr = datestr(t,'HH:MM:SS');
        
            if Device == "LAB"
           filenameLAB = sprintf("C:\\Users\\Quantum Dot\\Desktop\\Bera Yavuz - ANC300 Movement and Images\\QD_History_Logs\\QD_History_Text_Files\\%s",HistoryTextFileName);
           elseif Device == "HOME_PC"
           filenameLAB = sprintf("C:\\Users\\yavub\\OneDrive\\Desktop\\QD_History_logs\\All_QD_History_Files\\%s",HistoryTextFileName);
           elseif Device == "MAC"
           filenameLAB = "/Users/bera_yavuz/Desktop/Testing dumb thing folder/Fast_Movement_Settings.txt";
           else
               fprintf("invalid Device name try again")
            end


             if Device == "LAB" 
            filename = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\XY_Factor_Images\XY_Factor_History.txt";
            elseif Device == "HOME_PC"
            filename = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\XY_Factor_History.txt"; 
            else
                fprintf("invalid Device name try again")
            end
    

            switch InputSource 

                case "Write" 
                fileID = fopen(filename, 'a+');
                % Check if the file opened successfully
                if fileID == -1
                    error('Failed to open the file.');
                end

                % Get the current timestamp
                currentDateTime = datetime('now');
                currentDateTimeText = string(currentDateTime);
               
                % Write the timestamp and information to the file
                fprintf(fileID, '%s\n================================================\nCurrent Fast Movement X Step: %s\nnCurrent Fast Movement Y Step: %s\n\n\n', currentDateTimeText,XY_factor);

                
                case "Read"
                fileID = fopen(filename, 'rt');
                if fileID == -1
                    error('Failed to open the file.');
                end
                
                % Read the entire file into a single string
                fileContent = fread(fileID, '*char')';
                
                % Split the content into lines
                lines = strsplit(fileContent, newline);
                
                % Initialize variables to store the last three non-empty lines
                lastLine = '';
                secondLastLine = '';
               
                
                nonEmptyCount = 0;  % Counter for non-empty lines
                
                % Loop through the lines in reverse order
                for i = length(lines):-1:1
                    if ~isempty(strtrim(lines{i}))
                        nonEmptyCount = nonEmptyCount + 1;
                        if nonEmptyCount == 1
                            lastLine = lines{i};
                        elseif nonEmptyCount == 2
                            secondLastLine = lines{i};
                              % Stop after finding the second-to-last non-empty line
                        end
                    end
                    % Split the string into words
                    words = strsplit(lastLine);
                    % Extract the last word
                    Y_Factor = words{end};
                    % Split the string into words
                    words = strsplit(secondLastLine);
                    % Extract the last word
                    X_Factor = words{end};
                end
                
                % Convert the last and third-to-last lines to numbers
                Read_XY_factor = struct("X_factor",str2num(X_Factor),"Y_factor",str2num(Y_Factor)); 
                   
            end
            if Read_Write == "Write"
                Read_XY_factor = ''; 
            end
            fclose(fileID);
        
        end
      
        % Useful Convenience Functions  
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        function AddPathFunc(obj,Device)
             % adds needed paths to device 
            % Inputs:
                % Device - by specificying device user is able to add all paths needed at once (MAC,HOME_PC,LAB) 
            % Outputs:
                % Function has no variable that needs to be return 
            
            switch Device
                case "MAC"
                    directoryPath_AxesImg = "/Users/bera_yavuz/Desktop/GitHub ANC300 Project/Axes Images";
                    directoryPath_GoodImg = "/Users/bera_yavuz/Desktop/GitHub ANC300 Project/Photos Good For Testing";
                    directoryPath_Funcs = "/Users/bera_yavuz/Desktop/GitHub ANC300 Project";
                    direc_imgs = "/Users/bera_yavuz/Desktop/Quantum Dot Project Github/Bera-Yavuz-Quantum-Dot-Project-/Scripts For Testing_Debugging/Testing Images";
                    addpath(directoryPath_AxesImg,directoryPath_GoodImg,directoryPath_Funcs,direc_imgs);

                case "LAB"
                    directoryPath_Funcs = "C:\Users\Quantum Dot\Desktop\Bera_Yavuz_GitHub\AttoCube-Project-Stuff";
                    directoryPath_LED = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\LED_Find_Images"; 
                    directoryPath_Scripts = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools";
                    addpath(directoryPath_Funcs,directoryPath_LED,directoryPath_Scripts);
              
                case "HOME_PC"
                    directoryPath_AxesImg = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\Axes Images";
                    directoryPath_GoodImg = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff\Photos Good For Testing";
                    directoryPath_Funcs = "C:\Users\yavub\OneDrive\Desktop\ANC300 Project GitHub\AttoCube-Project-Stuff"; 
                    addpath(directoryPath_AxesImg,directoryPath_GoodImg,directoryPath_Funcs);
            end
        end

        function Current_QD_ID = Specific_QD_Movement(obj,ANC300,InquiredQD,SpectrometerPlot,pyueye_initialization_return)
            % Parameters
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            % X and y factors
            Read_XY_factor = XY_Factor_Identifier(obj,"","Read","LAB");
            x_factor = Read_XY_factor.X_factor; %(Steps/pixels units)
            y_factor = Read_XY_factor.Y_factor; %(Steps/pixels units)

            Frequency = 20; 
            % Adding needed pathways depending on Device - for lab purposes always use device name LAB 
            AddPathFunc(obj,"LAB")
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            % Initial Movement Towards Start
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            PhotoType = "SteppingPhoto";
            readQD_position = QD_tracking_N_Identification(obj,"","","Read","LAB","");
            QD_counter = readQD_position.lastLine; 
            CurrentQD =  QD_counter;
            [StartingQD,StartingQD_rotated,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = Precision_Locking(obj,ANC300,PhotoType,QD_counter,pyueye_initialization_return);
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            % calculating number of rows and column alongside direction needed to travel to inquired QD 
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            QD_loca_diff = CurrentQD - InquiredQD; 

            if QD_loca_diff(1) > 0 & QD_loca_diff(2) < 0

            % direction for each respectiive axis 
            direction_row = "topright";
            direction_column = "bottomright"; 

            % Intervals for each respective axis 

            row_interval = [-1 0];
            column_interval = [0 1];

            elseif QD_loca_diff(1)  > 0 & QD_loca_diff(2) > 0 

            % direction for each respectiive axis 
            direction_row = "topright";
            direction_column = "topleft";

            % Intervals for each respective axis 
            row_interval = [-1 0];
            column_interval = [0 -1];

            elseif QD_loca_diff(1)  < 0 & QD_loca_diff(2) < 0 

            % direction for each respectiive axis 
            direction_row = "bottomleft";
            direction_column = "bottomright"; 

            % Intervals for each respective axis 
            row_interval = [1 0];
            column_interval = [0 1];

            elseif QD_loca_diff(1)  < 0 & QD_loca_diff(2) > 0 
            
            % direction for each respectiive axis     
            direction_row = "bottomleft";
            direction_column = "topleft";

            % Intervals for each respective axis 
            row_interval = [1 0];
            column_interval = [0 -1];

            elseif QD_loca_diff(1) == 0 | QD_loca_diff(2) == 0
                if QD_loca_diff(1) < 0 
                direction_row = "bottomleft";
                row_interval = [1 0];
                elseif QD_loca_diff(1) > 0 
                direction_row = "topright";
                row_interval = [-1 0];
                elseif QD_loca_diff(2) < 0 
                direction_column = "bottomright"; 
                column_interval = [0 1];
                elseif QD_loca_diff(2) > 0 
                direction_column = "topleft"; 
                column_interval = [0 -1];
                else
                    fprintf("No movement")
                end

            end

            % number of rows and columns to move 
            num_rows = abs(QD_loca_diff(1));
            num_columns = abs(QD_loca_diff(2));
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


            %  Movement Towards QD  
            % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            PhotoType = "SteppingPhoto"; 

            for column = 1: num_columns
                % find next point along the column 
                [~,~,~,DistanceBetweenPoints_Rotated] = OptimizedRasterScan(obj,StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_column);

                % Moving to the next QD along column
                fprintf("QD Coord Difference: X = %.3f | Y = %.3f\n",DistanceBetweenPoints_Rotated)
                Dual_ANC300_Movement(obj,DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_column,ANC300,Frequency,x_factor,y_factor)

                % Making sure LED is exactly on the dot 
                QD_counter = QD_counter + column_interval; 
                QD_tracking_N_Identification(obj,QD_counter,"Auto","Write","LAB",direction_column); 
                [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = Precision_Locking(obj,ANC300,PhotoType,QD_counter,pyueye_initialization_return); % trying to land on the exact dot 
                fprintf("Currently on [%d %d] QD\n",QD_counter);
                Current_QD_ID = QD_counter; 
            end

            for row = 1:num_rows
                % finding next point along the row 
                [~,~,~,DistanceBetweenPoints_Rotated] = OptimizedRasterScan(obj,StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_row);
            
                % Moving to the next QD along row 
                fprintf("QD Coord Difference: X = %.3f | Y = %.3f\n",DistanceBetweenPoints_Rotated)
                Dual_ANC300_Movement(obj,DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_row,ANC300,Frequency,x_factor,y_factor)

                % Making sure LED is exactly on the dot 
                QD_counter = QD_counter + row_interval; 
                QD_tracking_N_Identification(obj,QD_counter,"Auto","Write","LAB",direction_row); 
                [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = Precision_Locking(obj,ANC300,PhotoType,QD_counter,pyueye_initialization_return); % trying to land on the exact dot 
                fprintf("Currently on [%d %d] QD\n",QD_counter);
                Current_QD_ID = QD_counter; 
            end

            fprintf("User has arrived to wanted QD: [%d %d]\n",InquiredQD)

            if SpectrometerPlot == "No Plot"
                fprintf("No plot was created"); 
            elseif SpectromerPlot == "Make Plot"
                % Spectroscope Image and Plot
                pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d]',QD_counter); 
                py.asi_func.snap_image(asi_initialization_return, pyrun_file_text_Spec)
            else
                error("Invalid option for SpectrometerPlot please try again")
            end


            % Exit 
            py.pyueye_func.exit_camera(pyueye_initialization_return)
            % ------------------------------------------------------------------------
        end

        function date_str = Fetch_Date(obj)
        %FetchDate Creates date string with leading zeros in DD_MM_YYYY format
        %   Returns date_str: String with zero-padded date (e.g., "20_01_2025")
        %
        %   Example:
        %   date_str = getFormattedDate();  % Returns "20_01_2025" for Jan 20, 2025
        %
        %   Uses current date by default. Modify datetime call for specific dates.
            % Create a datetime object with zero-padded format
            dt = datetime('now', 'Format', 'dd MM yyyy'); % For current date
            % dt = datetime(2025,1,20, 'Format', 'dd MM yyyy'); % For a specific date
            
            % Split the formatted date string into components
            date_str = string(dt);
            split_str = split(date_str, ' ');
            
            % Assign to variables with leading zeros
            day_str = split_str(1);    % Zero-padded day (e.g., "05")
            month_str = split_str(2);  % Zero-padded month (e.g., "02")
            year_str = split_str(3);   % Full year (e.g., "2025")
            
            % putting all the strings together for the date
            date_str = sprintf("%s_%s_%s",day_str,month_str,year_str); 
        end

        % Camera Related Functions (setting up, snapping, etc) 
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        function [vid_ASI,src_ASI,vid_UI,src_UI,camInfo] = ASI_UI_CameraInit(obj,ASI_Settings,UI_Settings)
        %ASI_UI_CAMERAINIT Initialize and configure ASI and UI cameras for imaging
        %   This function detects connected cameras, identifies ASI and UI devices,
        %   establishes video connections, and configures camera parameters.
        %
        %   Inputs:
        %   obj - Parent object reference (for class-based implementation)
        %   ASI_Settings - Struct containing configuration parameters for ASI camera:
        %       desired_exposure_time: Exposure time in seconds (converted to log2 value)
        %       Exposure: Gain value (0-100)
        %       Gamma: Gamma correction value
        %       Brightness: Brightness value (0-100)
        %   UI_Settings - Reserved for future UI camera configuration (currently unused)
        %
        %   Outputs:
        %   vid_ASI - Video input object for ASI camera
        %   src_ASI - Source object for ASI camera properties
        %   vid_UI - Video input object for UI camera
        %   src_UI - Source object for UI camera properties
        %
        %   Example:
        %   ASI_Config = struct('desired_exposure_time', 0.1, 'Exposure', 50, ...
        %                      'Gamma', 1.0, 'Brightness', 75);
        %   UI_Config = struct();
        %   [vidASI, srcASI, vidUI, srcUI] = ASI_UI_CameraInit([], ASI_Config, UI_Config); 

        % finding current imaging devices connected 
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
        vid_ASI = videoinput(adaptor, ASI_Device_ID,camInfo.DeviceInfo(ASI_Device_ID).DefaultFormat); % function can take a third input to specify formatting (ASK Sreesh) 
        src_ASI = getselectedsource(vid_ASI);
        %all_props_ASI = propinfo(vid_ASI); shows all settings user can change 
        
        % establishing connection and parameters of UI Device 
        vid_UI = videoinput(adaptor, UI_Device_ID);
        src_UI = getselectedsource(vid_UI);
        % all_props_UI = propinfo(vid_UI); % shows all settings user can change 
        
        % Setting correct Parameters for ASI camera
        if numel(fieldnames(ASI_Settings)) < 4 %checking that enough settings exist for everything to be properly assigned 
            fprintf("Not enough settings have been assigned to the ASI_Settings Struct")
            return
        end

        proper_exposure_setting = round(log2(ASI_Settings.desired_exposure_time)); %2^-15 seconds to 2^11 seconds 
        set(src_ASI,"ExposureMode","manual", "Exposure", proper_exposure_setting); %exposure is calculated as 2^n seconds by the camera (n = [-15 11])
        set(src_ASI,"GainMode", "manual" ,"Gain", ASI_Settings.Gain); % setting gain 
        set(src_ASI, "GammaMode","manual","Gamma",ASI_Settings.Gamma); % setting gamma
        set(src_ASI,"Brightness",ASI_Settings.Brightness); % setting brightness
        
        end
        
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
      
                    % Capture and save background images
                    start(vid_ASI)
                    wait(vid_ASI) % wait for frames to be acquired
                    background_img = getdata(vid_ASI,vid_ASI.FramesPerTrigger);

                    % get average of background and store it 
                    background_img = mean(background_img, 4);
                    if size(background_img,3) == 3 % checks if image is rgb
                        background_img = rgb2gray(background_img);
                    end
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

                    
                    % Snapping a photo and gray scaling it 
                    start(vid_ASI)
                    wait(vid_ASI) % wait for frames to be acquired
                    Emission_Reading_Img = getdata(vid_ASI,vid_ASI.FramesPerTrigger);

                    Emission_Reading_Img = mean(Emission_Reading_Img, 4);
                    if size(Emission_Reading_Img,3) == 3 % checks if image is rgb
                        Emission_Reading_Img = rgb2gray(Emission_Reading_Img);
                    end
                    % grab background average
                    background_Img = data.background_img; 
      
                    
                    % Auto-size vertical window
                    [height,width] = size(Emission_Reading_Img); 
                    central_row = height/2; 
                    window_size = round(height*0.15); % 40% below and above the central_row 
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
                            if ~exist(fullfile(latestFolder,"Raw_Imgs"),'dir')
                                mkdir(fullfile(latestFolder,"Raw_Imgs"))
                            end
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
                     fprintf(file,'%s \t %s\n', 'Wavelen. (nm)','Count Average');
                    fprintf(file,'%.4f \t %u\n',data_matrix);

            end
        end
    
        function [UI_Position_Img] = UI_Snap_Img(obj,vid_UI,src_UI,SaveImg,QD_ID) 
        % UI_Snap_Img - Captures and saves an image from the UI camera.
        %
        % Syntax:
        %   [UI_Position_Img] = UI_Snap_Img(obj, vid_UI, src_UI, SaveImg, QD_ID)
        %
        % Description:
        %   This function captures an image from the UI camera, which is used for position tracking. 
        %   The captured image is optionally saved based on the `SaveImg` parameter.
        %
        % Inputs:
        %   obj       - Object instance (not explicitly used in this function but required for class methods)
        %   vid_UI    - Video input object for the UI camera
        %   src_UI    - Source properties of the UI camera (not used in this function)
        %   SaveImg   - String flag ("Yes" or "No") indicating whether to save the captured image
        %   QD_ID     - Identifier for the Quantum Dot sample or NanoWire chip position
        %
        % Outputs:
        %   UI_Position_Img - Captured image from the UI camera
        %
        % Functionality:
        %   - Captures an image using the UI camera.
        %   - Saves the image to a predefined directory if `SaveImg` is "Yes".
        %   - The saved image is named based on the `QD_ID` and timestamp.
        %
            % turning off vertical and horizontal flip 

            
            % fetching today's date
            date = Fetch_Date(obj); 
            date_test =  date + "_Test"; 
            
            % Defining the pathway for where the UI Imgs end up 
            pathway_main = "c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\" + date_test +"\Position_uEYE"; 
            

            % Capturing a snap of the UI camera
            UI_Position_Img = getsnapshot(vid_UI);

            % Evaluating boolean to see if photo gets saved
            if SaveImg == "Yes"
                
                file_name = sprintf("[%d %d]_NanoWireChip_Pos.jpg",QD_ID);
                full_pathway = fullfile(pathway_main,file_name);
                imwrite(UI_Position_Img,full_pathway)

            elseif contains(SaveImg,"No")
            else

                file_name_debug = sprintf("%s.jpg",SaveImg); 
                full_pathway_debug = fullfile(pathway_main,file_name_debug); 
                imwrite(UI_Position_Img,full_pathway_debug)

            end

        end
    
        function [spectrum_sum] = ASI_Live_Feed_Snapping(obj,Emission_Reading_Img,bckgrnd_avg)
        
            
            % gray scaling photo
            if size(Emission_Reading_Img,3) == 3 % checks if image is rgb
                Emission_Reading_Img = rgb2gray(Emission_Reading_Img);
            end

            % Auto-size vertical window
            [height,width] = size(Emission_Reading_Img); 
            central_row = height/2; 
            window_size = round(height*0.4); % 40% below and above the central_row 
            valid_rows = max(1, central_row - window_size):min(height, central_row + window_size);

            % Sum spectrum with automated window
            spectrum_sum = sum(double(Emission_Reading_Img(valid_rows, :)), 1);
            
            spectrum_sum = spectrum_sum - bckgrnd_avg ;
        end
    
        % FSS Related functions
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        function stochastic_movement(obj, dir_x, dir_y, ANC300, Frequency)
            
            % Calculate the random perturbation (random step size)
            step_size_x = 4; % Random step size between 3 and 6
            step_size_y = 4; % Random step size between 3 and 6
        
        
            % step direction command 
            if dir_x == 1
            stochastic_dir_x = sprintf("stepu 1 %d",step_size_x); % Move up in x axis 
            else
            stochastic_dir_x = sprintf("stepd 1 %d",step_size_x); % Move down in x axis 
            end
            
            if dir_y == 1
            stochastic_dir_y = sprintf("stepu 2 %d",step_size_y); % Move up in y axis 
            else
            stochastic_dir_y = sprintf("stepd 2 %d",step_size_y); % Move down in y axis 
            end
        
            % sending movement command to piezo 
            fprintf(ANC300,stochastic_dir_x);
            StepQueue(obj,step_size_x,Frequency);
        
            fprintf(ANC300,stochastic_dir_y);
            StepQueue(obj,step_size_y,Frequency);
            
    
        end

        function Current_angle = FSS_Process(obj,ell_motor,QD_ID,vid_ASI,src_ASI,Spectrometer_Gratting,SaveRawImg,total_ange)

            % defining angle of rotation 
            angle = 0:5:total_ange;

            % for loop for rotating motor 
            for rot_count = 1:length(angle)

            % Defining the movement serial code for the rotation 
            angle_hxd = dec2hex(floor(mod(angle(rot_count),360)*39822/100), 8);
            input_str = "1ma" + angle_hxd; % "2" before ma is to be get from the ELLO software from thorlabs.
            
            % Commiting Command for movement 
            fprintf(ell_motor, input_str);
            
            
            Current_angle = sprintf('Angle: %d \n', angle(rot_count));
            fprintf(Current_angle)
            pause(1)
            
            % taking a snap at every respective angle 
            file_name = sprintf('FSS %d',angle(rot_count)); 
            ASI_Snap_Img(obj,vid_ASI,src_ASI,"Spectrometer",SaveRawImg,Spectrometer_Gratting,QD_ID,file_name);            
            end

        end
      
        function FSS_Folder_Creation(obj,QD)
            % Fetching today's date
            date = Fetch_Date(obj); % fetching today's date
            date_test =  date + "_Test"; 

    

            
            % Base folder path
            ASI_plots_directory = strcat(date_test, "\FSS_Measurements\Spectrometer_Plots_FSS\");
            text_file_pathway = fullfile("c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data", ASI_plots_directory);
            
            copy_number = 1;
            
            while true
                % Generate folder name with the current copy number
                folder_name = sprintf("QD_[%d %d]\\Set_%d", QD(1), QD(2), copy_number);
                folder_name = fullfile(text_file_pathway,folder_name); 
                
                % Check if the folder exists
                if ~exist(folder_name, 'dir')
                    mkdir(folder_name); % Create the folder if it does not exist
                    break; % Exit the loop once a new folder is created
                else
                    copy_number = copy_number + 1; % Increment the copy number
                end
            end
            
            % Final assigned folder path
            Specific_QD_FSS = folder_name;
            disp("Created folder: " + Specific_QD_FSS);
        end
        
        function latestFolder = find_latestFolder(obj,QD)
            % Fetching today's date
            date = Fetch_Date(obj); % fetching today's date
            date_test =  date + "_Test"; 
            

            
            % Base folder path
            ASI_plots_directory = strcat(date_test, "\FSS_Measurements\Spectrometer_Plots_FSS");
            Specific_QD_dir = sprintf("QD_[%d %d]",QD);
            text_file_pathway = fullfile("c:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data", ASI_plots_directory,Specific_QD_dir);

            parentDir = text_file_pathway;
            folders = dir(parentDir);
                folders = folders([folders.isdir]); % Keep only directories
                folders = folders(~ismember({folders.name}, {'.', '..'})); % Remove '.' and '..'
                
                if isempty(folders)
                    error('No folders found in the specified directory.');
                end
            
                % Sort by creation date (latest first)
                [~, idx] = sort([folders.datenum], 'descend');
                latestFolder = fullfile(parentDir, folders(idx(1)).name);
                
        end
    
    end
end