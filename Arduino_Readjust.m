% Arduino Angle Correction Code
Arduino = Arduino_Init; 
angle = read_encoder(Arduino); 


function Arduino = Arduino_Init()
    %connect to Arduino
    %serialPort = "COM8";
    % try
    %     Arduino = serialport(serialPort,'TimeOut',1,'BaudRate',9600);% Opening serial port to arduino
    % 
    serialPort = 'COM8';
    Arduino = serial(serialPort,'TimeOut',1,'BaudRate',9600);
    try
        fopen(Arduino);      % Opening serial port to arduino
        %handles.serialConnect = serialConnect;
    catch exception
        exception; %#ok<*VUNUS>
        exception.message;
        errordlg('Motor not connected','Connection Error')
        
    end
end


function angle = read_encoder(Arduino)
% read angle from encoder 
    flushinput(Arduino); % flush previously buffered info 
    try
        fprintf(Arduino,'%s','angle');
    catch
        error('Problem connecting to Arduino Board.');
    end
    pause(0.25);
    angle = str2num(fscanf(Arduino,'%s')); % read the current angle of the Arduino
end

% call this function to turn the grating
function sendArduino(mode,dir,size,Arduino)
    
    ArduinoFeed = [dir,mode,size];  % steps in "size" in the scan direction
    fprintf(Arduino,join(ArduinoFeed));
end


% No Idea what is happening in this code
function calibrateGratingSwitch(gratingPosition, newGratingPosition, originalPeak, originalWavelength, handles, hObject)

    centerWavelength = handles.centerWavelength;

    switch gratingPosition
        case 1
            wavelengthList = handles.wavelengths1;
        case 2
            wavelengthList = handles.wavelengths2;
        case 3  
            wavelengthList = handles.wavelengths3;
    end
    
    switch newGratingPosition
        case 1
            newWavelengthList = handles.wavelengths1;
        case 2
            newWavelengthList = handles.wavelengths2;
        case 3  
            newWavelengthList = handles.wavelengths3;
    end
    

    [~,targetIndex] = min(abs(centerWavelength-wavelengthList));
    [~,newTargetIndex] = min(abs(centerWavelength-newWavelengthList));

    orientationList = handles.orientation;
    difference = orientationList(newTargetIndex) - orientationList(targetIndex);

    fprintf('Difference: ');
    fprintf(num2str(difference));
    fprintf('\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % convert sign of difference to tilting direction
    if difference < 0 % the target ori. is at the cw of current ori.
        difference = difference + 1;
        dir = 'r'; % turn cw
    elseif difference > 0 % the target ori. is at the ccw of current ori.
        difference = difference - 1;
        dir = 'l'; % turn ccw
    end
    steps = abs(difference*3200); % this gives the steps required in sixteenth microstep mode, 3200 in sixteenth microstep mode = 1 deg of grating tilt
    % sorting the distance in angle into appropriate modes and steps
    full_step  = floor(steps/16.0);       % 16 sixteenth steps = 1 full step
    steps = mod(steps,16.0);               % update diff to be unresolved by full steps
    half_step = floor(steps/8.0);        % 8 sixteenth steps = 1 half step
    steps = mod(steps,8.0);               % update diff to be unresolved by half steps
    quarter_step = floor(steps/4.0);     % 4 sixteenth steps = 1 quarter step
    steps = mod(steps,4.0);               % update diff to be unresolved by quarter steps
    eighth_step = floor(steps/2.0);      % 2 sixteenth steps = 1 eighth step
    steps = mod(steps,2.0);               % update diff to be unresolved by eighth steps
    sixteenth_step = ceil(steps);        % the remainder of the distance will be covered by sixteenth steps with rounding up

    if full_step ~= 0.0
        sendArduino('f',dir,num2str(full_step),handles,hObject);
    end
    if half_step ~= 0.0
        sendArduino('h',dir,num2str(half_step),handles,hObject);
    end
    if quarter_step ~= 0.0
        sendArduino('q',dir,num2str(quarter_step),handles,hObject);
    end
    if eighth_step ~= 0.0
        sendArduino('e',dir,num2str(eighth_step),handles,hObject);
    end
    if sixteenth_step ~= 0.0
        sendArduino('s',dir,num2str(sixteenth_step),handles,hObject);
    end

    handles.tilted = 1;
    guidata(hObject,handles);
    
    fprintf(num2str(0.45*abs(difference)));
    fprintf('\n');
    pause(0.45*abs(difference));
    
    [peak, ~] = fetchPeak(hObject);
    i=0;
    while peak < 0.7*originalPeak % Repeat until peak is detected in camera
        sendArduino('f',dir,'25',handles,hObject); %Turn 25 steps
        pause(0.5);
        [peak, wavelength] = fetchPeak(hObject); % Read Peak
        i=i+1;
        if i > 15
            error('The grating did not switch. Please restart.');
            break
        end
    end
    
    [~,originalTargetIndex] = min(abs(originalWavelength-newWavelengthList));
    [~,newTargetIndex] = min(abs(wavelength-newWavelengthList));
    
    fprintf('Original wavelength: ');
    fprintf(num2str(newWavelengthList(originalTargetIndex)));
    fprintf('\n');
    
    fprintf('wavelength: ');
    fprintf(num2str(newWavelengthList(newTargetIndex)));
    fprintf('\n');
    
    difference = orientationList(newTargetIndex) - orientationList(originalTargetIndex); % Re-position to correct wavelength
    
    fprintf('Difference: ');
    fprintf(num2str(difference));
    fprintf('\n');
    
    tiltGrating(difference, handles, hObject);
end