% Define constants
c = 299792458; % speed of light in m/s
h = 6.62607015e-34; % plank's constant
fit_points = 1e3;% number of points to perform fit analysis on
C = colororder; % get color order for plotting

% defining qd ID
qd = [169 5];
power = "3_13"; 

%Initialise variables
Del_Lambda = nan(3,2); % wavelength
Del_f = nan(3,2); % frequency
FSS_delta = nan(3,2); % fine structure splitting from zero
FSS_std_dev = nan(3,2); % fine structure splitting standard deviation

% search current file for all needed things
desired_folder = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\28_02_2025_Test\FSS_Measurements\Spectrometer_Plots_FSS\QD_[169 5]\Set_2"; % put in the name of the folder with the desired angles 

fileList = dir(fullfile(desired_folder, '*.txt*')); % Find all .txt files in latestFolder
fileList = fileList(~[fileList.isdir]); % Remove directories from list
[~, idx] = sort([fileList.datenum]); % Sort files by date
fileList = fileList(idx); % Organize file list
numberOfFiles = length(fileList); % Define variable for total number of files

% Intiialse matrix for X, XX and HWP
j=1;

X_peaks_locs = nan(length(numberOfFiles),1);
XX_peaks_locs = nan(length(numberOfFiles),1);
hwp_angle_deg = nan(length(numberOfFiles),1);
end_angle = 5 * (numberOfFiles-1); 
%For angle not defined at start of filename define hwp outside for loop
hwp_angle_deg = 0:5:end_angle; % define HWP rotation angle in degrees
hwp_angle_deg = reshape(hwp_angle_deg,[length(hwp_angle_deg),1]);% reshape to column

N_vals = 100; % defining the number of initial points to average for more clear data
X_ID = 2; % assign which ever peak is the exciton 
XX_ID = 1; % assign which ever peak is the biexciton
        
 

for i = 1:numberOfFiles
    fileName = fileList(i).name;% go to i-th file
    data = readtable(fullfile(desired_folder,fileName),'HeaderLines', 7, 'PreserveVariableNames', true);% read data from txt. file, skip fist 7 lines
    % data = table2array(data(:,1:end-1)); % convert data from table to array, will exclude the last colunm with -1 (sometimes an empty column appears and causes an issue)
    data = table2array(data(:,1:end)); % convert data from table to array, will exclude the last colunm with -1 (sometimes an empty column appears and causes an issue)
    lambda = data(:,1);% extract wavelength from 1st column
    qd_emission = mean(data(:,2:end),2); % extract counts and calculte mean from no. figures avergaed over
    bckgrnd_val = sum(qd_emission(1:N_vals))/N_vals; 
    debug_file = i; % if code fails can see which file it failed on
      if i==1
        [pks, locs] = findpeaks(qd_emission, lambda,'MinPeakHeight',80,'NPeaks',4,'SortStr','descend','MinPeakDistance',0.2); % find peaks from data using find peak function uses y then x
        [~, idx] = sort(locs); % sort order of peak in order of wavelength, will default to highest no. counts if not.
        locs = locs(idx);% create variable for peak location (lambda)
        pks = pks(idx);% sort peak according to location (lambda)
        X_loc = locs(X_ID);% define first peak location as X
        XX_loc = locs(XX_ID);% define second peak location as XX
        lambda_X_initial = X_loc; % define the initial varaible of lambda to be location of X
        % Peak is not a smooth curve so define range (50pm) from which to apply fit for each peak
        X_low = X_loc-0.050; % X min range value
        X_high = X_loc+0.050; % X max range value
        XX_low = XX_loc-0.050; % as above for XX
        XX_high = XX_loc+0.050;
        % plot figure to verify fitting range is good.
        figure(1); clf;
            plot(lambda, qd_emission); grid on; hold on;
            plot([X_low X_low], [0 pks(X_ID)],'r');
            plot([X_high X_high], [0 pks(X_ID)],'r');
            plot([XX_low XX_low], [0 pks(XX_ID)],'g');
            plot([XX_high XX_high], [0 pks(XX_ID)],'g');
    end

    xdata_X = lambda(lambda>X_low & lambda<X_high); % new varibale to define the peak range (wavelength)
    ydata_X = qd_emission(lambda>X_low & lambda<X_high); % new varibale to define the peak range (counts)
    [x_fit, y_fit, fit_model, gof] = fit_lorentzian(xdata_X, ydata_X, fit_points, bckgrnd_val);% custom function written by Rubayet   , gof = goodness of fit 
    %[x_fit, y_fit, fit_model, gof] = fit_gaussian(xdata_X, ydata_X, fit_points, bckgrnd_val);% custom function written by Rubayet   , gof = goodness of fit 
    [pks, locs] = findpeaks(y_fit, x_fit,'MinPeakHeight',bckgrnd_val);% find peak of fitted function
    X_peaks_locs(i,:) = locs; % find location of max point

    xdata_XX = lambda(lambda>XX_low & lambda<XX_high);% as above but for XX
    ydata_XX = qd_emission(lambda>XX_low & lambda<XX_high);
    [x_fit2, y_fit2, fit_model2, gof2] = fit_lorentzian(xdata_XX, ydata_XX, fit_points, bckgrnd_val); 
    %[x_fit2, y_fit2, fit_model2, gof2] = fit_gaussian(xdata_XX, ydata_XX, fit_points, bckgrnd_val); 
    [pks, locs] = findpeaks(y_fit2, x_fit2,'MinPeakHeight', bckgrnd_val);
    XX_peaks_locs(i,:) = locs;

end

% save .mat file with all the extracted parameters from the plots related
% to the X and XX.
% Generate dynamic .mat file name
matFileName = sprintf('FSS_%d_%d_QD_analysis_P%s_uW_V3.mat', qd, power);
matFileName = fullfile(desired_folder,matFileName); 

% Save variables to the dynamic .mat file
save(matFileName, 'xdata_X', 'ydata_X', 'X_peaks_locs', 'xdata_XX', 'ydata_XX', 'XX_peaks_locs', 'hwp_angle_deg', 'X_loc', 'lambda_X_initial');
%%
%convert peak location from lambda (nm) to energy (meV)
E_X = h.*c./(X_peaks_locs.*1e-9)./(1e-3*1.602e-19); % in meV
E_XX = h.*c./(XX_peaks_locs.*1e-9)./(1e-3*1.602e-19);  % in meV

% X - XX energies
FSS = 1e3.*(E_X-E_XX);  % in mueV
FSS = FSS-mean(FSS);   % in mueV, to oscillate around zero value for plotting
hwp_angle_rad = hwp_angle_deg.*pi./180; % convert degrees to radians

%Plot X and XX energy as HWP is rotated
f4=figure(4);clf;box on;
set(gcf,'paperpositionmode','auto','color','w');
subplot(2,1,1)
plot(hwp_angle_deg,X_peaks_locs,'ob')
xlabel('HWP angle (degrees)')
ylabel('\lambda_X (nm)')
subplot(2,1,2)
plot(hwp_angle_deg,XX_peaks_locs,'or')
xlabel('HWP angle (degrees)')
ylabel('\lambda_{XX} (nm)')
% ylim([894.8965 894.8972])

xdata = hwp_angle_rad;% redefine x and y data
ydata = FSS;
[x_fit, y_fit, fit_model3, gof] = fit_sine(xdata, ydata, fit_points);% custom function from Rubayet
x_fit_deg = linspace(x_fit(1).*180./pi, x_fit(end).*180./pi, fit_points);% convert x data from rad to deg.

% generate error bar using matlab fit_model
ci = confint(fit_model3,0.95);
FSS_std = diff(ci(:,1))/2;

% generate plot
f2=figure(2); clf;%box on;
plot(2.*hwp_angle_deg,FSS-fit_model3.d,'o','MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:)); % plot data
hold on; grid on;
plot(2.*x_fit_deg, y_fit-fit_model3.d,"Color",C(3,:))% plot fit
xlabel('Polarization rotation angle (rad.)');
ylabel('$E_{X}-E_{XX} - \overline{(E_{X}-E_{XX})} (\mu eV)$','Interpreter','latex');

title(sprintf('QD (%d,%d), FSS = %.3f \\pm %.2f \\mueV', qd, fit_model3.a, FSS_std));
% del_lambda = 1e3.*(X_loc-lambda_X_initial); % (in pm), change due to gas deposition
% del_f = c.*(1/lambda_X_initial-1/X_loc);% convert change to frequency
% title("\Delta\lambda_X = "+ num2str(del_lambda, 3)+" pm, \Deltaf = "+num2str(del_f, 3)+" (GHz),  FSS = "+num2str(fit_model3.a, 3)+"\pm"+num2str(FSS_std, 2)+" \mueV")
ax = gca;
ax.XTick = 0:90:360;
ax.XTickLabel = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi', '5\pi/2', '3\pi', '7\pi/2', '4\pi'}; 
% ylim([-4 4])

% Del_Lambda(1,1) = del_lambda;
% Del_f(1,1) = del_f;
FSS_delta(1,1) = fit_model3.a;
FSS_std_dev(1,1) = FSS_std;

% save FSS_67_1_analysis.mat hwp_angle_deg FSS-fit_model3.d x_fit_Del_Lambda Del_f FSS_delta FSS_std_dev 

%% Plot nice figures
%Create handle for figure
f33=figure(33);clf;
fs = 12; %font size
pw = 12; %page width
ph = 8; %page height
%Configure Figure Size to single column width
set(f33,'unit','centimeters','color','w','position',[0 10 pw ph],...
'paperPositionMode','auto','ToolBar','none','MenuBar','none');
%Create Axis (note customised size)
ax=axes('units','centimeters','position',[1.2 1.2 pw-2 ph-2],'box','on',...
'linewidth',0.5,'fontsize',fs-2,'nextplot','add');
%Label Axis (xlbl+ylbl are handles to control position)
xlb=xlabel('HWP Rotation angle (rad.)','fontsize',fs,'interpreter','latex','units','centimeters');
pos=get(xlb,'position');set(xlb,'position',pos+[0.2 0.0 0]);
ax.XTick = 0:90:360;
ax.XTickLabel = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi', '5\pi/2', '3\pi', '7\pi/2', '4\pi'}; 
ylb=ylabel('$E_{X}-E_{XX} - \overline{(E_{X}-E_{XX})}~(\mu eV)$', 'fontsize', fs,'interpreter','latex','units','centimeters');
pos=get(ylb,'position');set(ylb,'position',pos+[0.2 0.0 0]);
ylim([-8 8])
plot(2.*hwp_angle_deg,FSS-fit_model3.d,'o','MarkerEdgeColor',C(1,:),'MarkerFaceColor',C(1,:)); % plot data
hold on; grid on;
plot(2.*x_fit_deg, y_fit-fit_model3.d,"Color",C(3,:))% plot fit
title(sprintf('QD (%d,%d), FSS = %.3f \\pm %.2f \\mueV', qd, fit_model3.a, FSS_std));

% Generate dynamic .png file name
pngFileName = sprintf('FSS_%d_%d_QD_P%s_uW_V3.png', qd, power);
pngFileName = fullfile(desired_folder,pngFileName); 

% Save the figure with the dynamic file name
print('-f33', '-dpng', '-r96', pngFileName);

disp(['Figure saved as: ', pngFileName]);