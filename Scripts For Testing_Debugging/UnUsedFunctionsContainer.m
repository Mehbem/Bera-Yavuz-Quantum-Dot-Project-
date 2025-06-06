classdef UnUsedFunctionsContainer
    methods
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
            angle = 44.8; % Pre Aug 16th 43.696474776653320; % use the angle found by the algorithm in the XYANC300Axes script 

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
        
        function UpdateText(obj,X_step_num,Y_step_num, time_to_pause,axis)

            switch axis 
                case "1"
                fprintf("---------------------------------------\nStepping Complete\n%d step(s) in X direction\nExpected time of Completion: %.3f seconds\n---------------------------------------\n",X_step_num,time_to_pause)
                case "2"
                fprintf("---------------------------------------\nStepping Complete\n%d step(s) in Y direction\nExpected time of Completion: %.3f seconds\n---------------------------------------\n",Y_step_num,time_to_pause)    
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
                if InputSource == "Auto"
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
     
    end
end

