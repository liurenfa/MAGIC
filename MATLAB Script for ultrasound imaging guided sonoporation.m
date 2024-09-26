clear all

% Frames_to_record = 1300;
% start_recording = 0;
markHandle = [];

delay_rec = 10; % delay in us between the Verasonics trigger out and the 
% start of the receive event. This value should be a bit shorter than the 
% travel time from the 250 kHz transducer to the fiber plus from the fiber 
% to the receiving transducer. 
% record_duration = 512; % in wavelength, multiply by 8 to get the number of 
% time samples recorded (i.e. 192 wavelengths gives 256*8 = 2048 time 
% samples). 

treat_pulse_length = 4000;
treat_pulse_num = 1;
% treat_pulse_duty = 0.1;

na_Bmode = 21;      % Set na = number of angles.
if (na_Bmode > 1), dtheta_Bmode = (20*pi/180)/(na_Bmode-1); else dtheta_Bmode=0; end % set dtheta to range over +/- 10 degrees.

na_CPS = 21;      % Set na = number of angles.
if (na_CPS > 1), dtheta_CPS = (20*pi/180)/(na_CPS-1); else dtheta_CPS=0; end % set dtheta to range over +/- 5 degrees.

na = max(na_Bmode,na_CPS);

PRF_Bmode = 50; % in Hz
PRF_CPS = 50; % in Hz

Nbframes = PRF_CPS*2; % 2s

%% Define system parameters.
Resource.Parameters.connector = 2;
Resource.Parameters.numTransmit = 128; 
Resource.Parameters.numRcvChannels = 128;  
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0; % 0 means no simulation, if hardware is present.
Resource.Parameters.startEvent = 1;
Resource.HIFU.externalHifuPwr = 1;
Resource.HIFU.extPwrComPortID = 'COM3';
% Resource.Parameters.fakeScanhead = 1;

%% Specify Trans structure array.
Trans.name = 'JRL9-4';
%Trans.name = 'CL15-7';
% Trans.name = 'JRL11-4';

% Trans.name = 'IP-104';
Trans.units = 'mm'; 
Trans = computeTrans(Trans);  
% Trans.frequency = 4.5;
Trans.maxHighVoltage = 60;  

P.startDepth = 5; % %50
P.endDepth = P.startDepth+50; %64

%% Specify PData structure array.

% TX(4).Apod = [zeros(1,14) ones(1,60-14) zeros(1,128-60)];

PData(1).PDelta = [0.5, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((127*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*127/2,0,P.startDepth]; % x,y,z of upper lft crnr.

% evalin('base','PData(1).Size(2) = ceil((max(find(TX(4).Apod))-min(find(TX(4).Apod)))/PData(1).PDelta(1));');
%         evalin('base','PData(1).Origin = [Trans.ElementPos(min(find(TX(4).Apod)),1),0,P.startDepth];');
%         evalin('base','PData(1).Region = computeRegions(PData(1));');

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048*na; %896*na*m*3; %64 1152
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = Nbframes;    
Resource.InterBuffer(1).numFrames = 1;   
Resource.ImageBuffer(1).numFrames = Nbframes;
Resource.ImageBuffer(2).numFrames = 20;
Resource.DisplayWindow(1).Title = 'L4-9 CPS';
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).pdelta = 0.15;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
% Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).splitPalette = 1;

% load customGamma_Dec2017
% Resource.DisplayWindow(1).Colormap = repmat(customGamma.Curve,1,3);
Resource.DisplayWindow(1).ReverseAxis = 0;

TPC(1).name = 'B-mode';
TPC(1).maxHighVoltage = 30;
TPC(1).hv = 2;

TPC(2).name = 'CPS';
TPC(2).maxHighVoltage = 50;
TPC(2).hv = 4;

TPC(5).name = 'Treat';
TPC(5).maxHighVoltage = 35;
TPC(5).hv = 1.6;
%% Specify Transmit waveform structure.  
TW(1).type = 'parametric';
TW(1).Parameters = [7.9,0.67,2,1];   % B-mode

TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency-1.8,0.67,2,1];   % positive polarity for CPS

TW(3).type = 'parametric';
TW(3).Parameters = [Trans.frequency-1.8,0.67,2,-1];   % negative polarity for CPS

TW(4).type = 'parametric';
TW(4).Parameters = [3.2,0.67,treat_pulse_length*2,1];   % positive polarity for CPS
%% Specify TX structure array.
vect_apod = ones(1,Resource.Parameters.numTransmit);
TX = repmat(struct('waveform', 1, ...
            'Origin', [0.0,0.0,0.0], ...
            'Apod', vect_apod, ... 
            'focus', 0.0, ...
            'Steer', [0.0,0.0], ...
            'Delay', zeros(1,Resource.Parameters.numTransmit)),1,na_Bmode+3*na_CPS+1); 

if fix(na_Bmode/2) == na_Bmode/2       % if na even
    P.startAngle_Bmode = (-(fix(na_Bmode/2) - 1) - 0.5)*dtheta_Bmode;
else
    P.startAngle_Bmode = -fix(na_Bmode/2)*dtheta_Bmode;
end
if fix(na_CPS/2) == na_CPS/2       % if na even
    P.startAngle_CPS = (-(fix(na_CPS/2) - 1) - 0.5)*dtheta_CPS;
else
    P.startAngle_CPS = -fix(na_CPS/2)*dtheta_CPS;
end

% B-mode
for n = 1:na_Bmode   % na transmit events
    angle = P.startAngle_Bmode + (n-1)*dtheta_Bmode;
    
    TX(n).Steer = [angle,0.0];
    TX(n).Delay = computeTXDelays(TX(n)); 
end

% CPS
for n = (1:3:3*na_CPS)+na_Bmode   % na transmit events
    angle = P.startAngle_CPS + (n-1-na_Bmode)/3*dtheta_CPS;
    
    TX(n).Steer = [angle,0.0];
    TX(n).Delay = computeTXDelays(TX(n)); 
    TX(n).Apod = zeros(1,Resource.Parameters.numTransmit);
    TX(n).Apod(1:2:127) = vect_apod(1:2:127);
    TX(n).waveform = 2;
    
    TX(n+1).Steer = TX(n).Steer ;
    TX(n+1).Delay = TX(n).Delay;
    TX(n+1).waveform = 3;
    
    TX(n+2).Steer = TX(n).Steer ;
    TX(n+2).Delay = TX(n).Delay; 
    TX(n+2).Apod = zeros(1,Resource.Parameters.numTransmit);
    TX(n+2).Apod(2:2:128) = vect_apod(2:2:128);
    TX(n+2).waveform = 2;
end

% TX(end).Apod = [zeros(1,14) ones(1,60-14) zeros(1,128-60)];
TX(end).Apod = ones(1,Resource.Parameters.numTransmit);
TX(end).focus = 0.0;
TX(end).waveform = 4;
%% Specify Receive structure arrays. 
InputFilterImg =[ -0.00113 +0.00000 -0.00116 +0.00000 +0.00549 +0.00000 +0.00720 ...
 +0.00000 -0.01419 +0.00000 -0.02640 +0.00000 +0.02606 +0.00000 ...
 +0.07816 +0.00000 -0.03671 +0.00000 -0.30786 +0.00000 +0.54108]; 

InputFilterPCD = [-0.00085 +0.00000 -0.00146 +0.00000 +0.00046 +0.00000 +0.00793 ...
 +0.00000 +0.01874 +0.00000 +0.02121 +0.00000 -0.00192 +0.00000 ...
 -0.06039 +0.00000 -0.14297 +0.00000 -0.21719 +0.00000 +0.75293];

delay_rec_lambda = delay_rec;
maxAcqLength = ceil(sqrt(P.endDepth^2 + (128*Trans.spacingMm*Trans.frequency*1e-3/Resource.Parameters.speedOfSound)^2)) ;
Receive = repmat(struct('Apod', ones(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'demodFrequency',2*(Trans.frequency-1.8),...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0),1, na*Resource.RcvBuffer(1).numFrames);%'InputFilter', InputFilterImg, ...

RcvProfile.PgaGain = 24; % 30
RcvProfile.LnaGain = 24; % 15, 18, 24
RcvProfile.LnaZinSel = 0; % 0 to 31
                    
% - Set event specific Receive attributes for each frame.
for i = 1:Nbframes
    for k = 1:na
        Receive((i-1)*na + k).framenum = i;
        Receive((i-1)*na + k).acqNum = k;
    end
end

Receive = [Receive Receive Receive];
for n = (1:2*na*Nbframes) + na*Nbframes
    Receive(n).mode = 1;
end
Receive(na*Resource.RcvBuffer(1).numFrames*3+1) = Receive(1);
            Receive(na*Resource.RcvBuffer(1).numFrames*3+1).startDepth = delay_rec_lambda;
            Receive(na*Resource.RcvBuffer(1).numFrames*3+1).endDepth = delay_rec_lambda+100;

%% Specify TGC Waveform structure.
TGC.CntrlPts = linspace(200,1023,8);%[0,141,275,404,510,603,702,782];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon structure arrays.
Recon(1) = struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:na_Bmode));
           
Recon(2) = struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:na_CPS)+na_Bmode);
Recon(3) = struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',na_CPS+na_Bmode+1);


% Define ReconInfo structures.               
ReconInfo = repmat(struct('mode','accumIQ', ...  %'replaceIQ' accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1,...
                   'scaleFactor',3), 1, na_Bmode+na_CPS+1);
% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ';
for j = 1:na_Bmode
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
end
ReconInfo(na_Bmode).mode ='accumIQ_replaceIntensity' ; %'replaceIntensity'

ReconInfo(na_Bmode+1).mode = 'replaceIQ';
for j = (1:na_CPS)+na_Bmode
        ReconInfo(j).txnum = (j-1-na_Bmode)*3+2+na_Bmode;
        ReconInfo(j).rcvnum = j-na_Bmode;
end
ReconInfo(na_Bmode+na_CPS).mode ='accumIQ_replaceIntensity' ; %'replaceIntensity'

ReconInfo(end).mode = 'replaceIntensity';
%% Specify Process structure array.
pers = 0;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',2.5,...            % pgain is image processing gain
                         'reject',0,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','log',...
                         'compressFactor',55,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
Process(2).classname = 'External';
Process(2).method = 'plot_range_image';
Process(2).Parameters = {'srcbuffer','none',...  % name of buffer to process.
                         'dstbuffer','image' ...
                         'dstbufnum',2, ...
                         'dstframenum',-2};
persf = 80;
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','signedColor',... % type of data to display.
                         'persistMethod','dynamic',...
                         'persistLevel',persf,...
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',0,...      % reject level
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',50,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

%% Specify SeqControl structure arrays.

% B-mode voltage
TPC1 = 1;
SeqControl(TPC1).command = 'setTPCProfile';
SeqControl(TPC1).argument = 1;
SeqControl(TPC1).condition = 'immediate';

% CPS voltage
TPC2 = 2;
SeqControl(TPC2).command = 'setTPCProfile';
SeqControl(TPC2).argument = 2;
SeqControl(TPC2).condition = 'immediate';

% Treat voltage
TPC3 = 3;
SeqControl(TPC3).command = 'setTPCProfile';
SeqControl(TPC3).argument = 5;
SeqControl(TPC3).condition = 'immediate';

TTNATX = 4;
SeqControl(TTNATX).command = 'timeToNextAcq';  
SeqControl(TTNATX).argument = 100;  % 

TTNFB = 5;
SeqControl(TTNFB).command = 'timeToNextAcq';  % PRF Bmode
SeqControl(TTNFB).argument = round(1e6/PRF_Bmode)-na_Bmode*SeqControl(TTNATX).argument;  

TTNFCPS = 6;
SeqControl(TTNFCPS).command = 'timeToNextAcq';  % PRF CPS
SeqControl(TTNFCPS).argument = round(1e6/PRF_CPS)-na_CPS*3*SeqControl(TTNATX).argument;  

RETMAT = 7;
SeqControl(RETMAT).command = 'returnToMatlab';

JUMPB = 8;
SeqControl(JUMPB).command = 'jump'; % jump back to B-mode
SeqControl(JUMPB).argument = 2;

JUMPCPS = 9;
SeqControl(JUMPCPS).command = 'jump'; % jump back to CPS
SeqControl(JUMPCPS).argument = 1;

JUMPTREAT = 10;
SeqControl(JUMPTREAT).command = 'jump'; % jump back to treatment
SeqControl(JUMPTREAT).argument = 1;

treat_noop = 11;
SeqControl(treat_noop).command = 'noop'; % jump back to treatment
SeqControl(treat_noop).argument = 5e3;%round(treat_pulse_length/TW(4).Parameters(1)*(1/treat_pulse_duty-1)/0.2);%treat_pulse_length/TW(4).Parameters(1)*(1/treat_pulse_duty-1) in us

TTNFCPSP = 12;
SeqControl(TTNFCPSP).command = 'timeToNextAcq';  % PRF CPS
SeqControl(TTNFCPSP).argument = SeqControl(TTNFCPS).argument*5;  
nsc = length(SeqControl)+1; % nsc is next count of SeqControl objects

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

%% Specify Event structure arrays.
n = 1;

Event(n).info = 'TPC(1)';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TPC1;
n = n+1;

%% B-mode
for i = 1:Nbframes
    for j = 1:na_Bmode     
        Event(n).info = 'B-mode';
        Event(n).tx = j;
        Event(n).rcv = (i-1)*na + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
    end
    Event(n-1).seqControl = [TTNFB,nsc];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)
        Event(n).seqControl = RETMAT;
    end
    n = n+1;
    
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        
Event(n).rcv = 0;      
Event(n).recon = 0;    
Event(n).process = 0; 
Event(n).seqControl = JUMPB;
n = n+1;

%% CPS mode 
n_CPS = n;

Event(n).info = 'TPC(2)';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TPC2;

SeqControl(JUMPCPS).argument = n;
n = n+1;
k = na*Nbframes;
for i = 1:Nbframes
    for j = 1:na_CPS
        Event(n).info = 'Amp=1/2';
        Event(n).tx = na_Bmode + (j-1)*3 + 1;
        Event(n).rcv = (i-1)*na + j   ;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
        
        Event(n).info = 'Amp=-1';
        Event(n).tx =  na_Bmode + (j-1)*3 + 2;
        Event(n).rcv = (i-1)*na + j + k;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
        
        Event(n).info = 'Amp=1/2';
        Event(n).tx =  na_Bmode + (j-1)*3 + 3;
        Event(n).rcv = (i-1)*na + j + 2*k;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
    end
    Event(n-1).seqControl = [TTNFCPS,nsc];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)
        Event(n).seqControl = RETMAT;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        
Event(n).rcv = 0;      
Event(n).recon = 0;    
Event(n).process = 0; 
Event(n).seqControl = JUMPCPS;
n = n+1;

%% CPS mode with treatment
n_treat = n;

Event(n).info = 'TPC(2)';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TPC2;


SeqControl(JUMPTREAT).argument = n;
n = n+1;
k = na*Nbframes;
for i = 1:Nbframes
    for j = 1:na_CPS
        Event(n).info = 'Amp=1/2';
        Event(n).tx = na_Bmode + (j-1)*3 + 1;
        Event(n).rcv = (i-1)*na + j   ;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
        
        Event(n).info = 'Amp=-1';
        Event(n).tx =  na_Bmode + (j-1)*3 + 2;
        Event(n).rcv = (i-1)*na + j + k;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
        
        Event(n).info = 'Amp=1/2';
        Event(n).tx =  na_Bmode + (j-1)*3 + 3;
        Event(n).rcv = (i-1)*na+ j + 2*k;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNATX;
        n = n+1;
    end
    Event(n-1).seqControl = [TTNFCPS,nsc];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)
        Event(n).seqControl = RETMAT;
    end
    n = n+1;   
    
    Event(n).info = 'plot';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;
    Event(n).info = 'display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1; 
end

    
   Event(n).info = 'TPC(3)';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = TPC3;
    n = n+1; 
    
for ii = 1:treat_pulse_num
    Event(n).info = 'noop';
    Event(n).tx =  0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = treat_noop;
    n = n+1;
    
    Event(n).info = 'Treat pulse';
    Event(n).tx =  na_Bmode+3*na_CPS+1;
    Event(n).rcv = na*Resource.RcvBuffer(1).numFrames*3+1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [TTNFCPS,nsc];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 3;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1; 
    
%     Event(n).info = 'plot';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 2;
%     Event(n).seqControl = 0;
%     n = n+1;
%     Event(n).info = 'display';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 3;
%     Event(n).seqControl = 0;
%     n = n+1; 
end
    
%     Event(n-1).seqControl = [TTNFCPSP,RETMAT];

Event(n).info = 'Jump back';
Event(n).tx = 0;        
Event(n).rcv = 0;      
Event(n).recon = 0;    
Event(n).process = 0; 
Event(n).seqControl = JUMPTREAT;
n = n+1;

wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end

frameRateFactor = 5;

%% User specified UI Control Elements
% Callback routine to start recording
% UI(1).Control =  {'UserA1','Style','VsToggleButton','Label','record'};
% UI(1).Callback = text2cell('%-UI#1Callback');
% import vsv.seq.uicontrol.VsSliderControl


UI(1).Control = {'UserB4','Style','VsButtonGroup','Title','Mode','NumButtons',2,'Labels',{'B-mode','CPS'}};
UI(1).Callback = text2cell('%-UI#1Callback');
% 
% UI(2).Control =  {'UserB1','Style','VsToggleButton','Label','foc. area'};
% UI(2).Callback = text2cell('%-UI#2Callback');

UI(3).Control =  {'UserC1','Style','VsToggleButton','Label','TREAT'};
UI(3).Callback = text2cell('%-UI#3Callback');

% - Left element Slider
UI(4).Control =  {'UserB3','Style','VsSlider','Label','Left range',...
        'SliderMinMaxVal',[Trans.ElementPos(1,1),Trans.ElementPos(Trans.numelements,1),Trans.ElementPos(1,1)],...
        'SliderStep',[0.5/(Trans.ElementPos(Trans.numelements,1)-Trans.ElementPos(1,1)),2.5/(Trans.ElementPos(Trans.numelements,1)-Trans.ElementPos(1,1))],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%LeftEleCallback');

% UI(4).Control = VsSliderControl('LocationCode','UserB3',...
%     'Label','Left element',...
%     'SliderMinMaxVal',[1,Trans.numelements,1+1],...
%     'SliderStep',[1/128 5/128],'ValueFormat','%3.0i',...
%     'Callback',@LeftEleCallback);

% - Right element Slider
UI(5).Control =  {'UserB2','Style','VsSlider','Label','Right range',...
        'SliderMinMaxVal',[Trans.ElementPos(1,1),Trans.ElementPos(Trans.numelements,1),Trans.ElementPos(Trans.numelements,1)],...
        'SliderStep',[0.5/(Trans.ElementPos(Trans.numelements,1)-Trans.ElementPos(1,1)),2.5/(Trans.ElementPos(Trans.numelements,1)-Trans.ElementPos(1,1))],'ValueFormat','%3.0f'};
UI(5).Callback = text2cell('%RightEleCallback');

% UI(5).Control = VsSliderControl('LocationCode','UserB2',...
%     'Label','Right element','SliderMinMaxVal',[1,Trans.numelements,Trans.numelements-1],...
%     'SliderStep',[1/128 5/128],'ValueFormat','%3.0i',...
%     'Callback', @RightEleCallback);


% Save all the structures to a .mat file.
save('Bmode_CPS_treat_L4_9_20231010');

return


%% **** Callback routines to be converted by text2cell function. ****
%-UI#1Callback - %%%%%%%%%%%%%%%%%%%%%%%% activation %%%%%%%%%%%%%%%%%%%%%%
% keyboard
Resource = evalin('base','Resource');
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
graycolormap = gray(256);
coppercolormap = copper(256);
switch UIState
   case 1  % B mode
        Control(1).Parameters = {'Parameters',1,'startEvent',1};
        Control(2).Command = 'set&Run';
        evalin('base',['Resource.Parameters.startEvent = ',num2str(1),';']);
        evalin('base','Resource.DisplayWindow(1).Colormap = gray(256);');
        Control(2).Parameters = {'DisplayWindow',1,'colormap',graycolormap};
%         Resource.DisplayWindow(1).Colormap = copper(256);
%         assignin('base','Resource',Resource);
        assignin('base','Control',Control);
   case 2  % CPS mode
        n_CPS = evalin('base','n_CPS');
        Control(1).Parameters = {'Parameters',1,'startEvent',n_CPS};
        Control(2).Command = 'set&Run';
        evalin('base',['Resource.Parameters.startEvent = ',num2str(n_CPS),';']);
        evalin('base','Resource.DisplayWindow(1).Colormap = copper(256);');
        Control(2).Parameters = {'DisplayWindow',1,'colormap',coppercolormap};
%         Resource.DisplayWindow(1).Colormap = copper(256);
%         assignin('base','Resource',Resource);
        assignin('base','Control',Control);
end

return
%-UI#1Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-UI#2Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% markHandle = evalin('base','markHandle');
% if isempty(markHandle)
%     bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
%     theta = linspace(0,2*pi,100)';
%     CT = 3.5*cos(theta);
%     ST = 3.6*sin(theta);
% end
% 
% if isempty(markHandle)
%     figure(bmodeFigHandle), 
%     hold on,
%     markHandle = plot(0+CT,10+ST,'--g');
%     hold off 
% end
% 
% switch UIState
%     case 0
%         markHandle.Visible = 'off';
%     case 1
%         markHandle.Visible = 'on';
% end
% assignin('base','markHandle',markHandle);
% return
%-UI#2Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-UI#3Callback - %%%%%%%%%%%%%%%%%%%%%%%% recording %%%%%%%%%%%%%%%%%%%%%%
Resource = evalin('base','Resource');
% TX = evalin('base','TX');
% Trans = evalin('base','Trans');
% P = evalin('base','P');
% PData = evalin('base','PData');
Process = evalin('base','Process');
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
graycolormap = gray(256);
copper0 = copper(256);
coppercolormap = grayscaleCPAmap;
coppercolormap(1:128,:) = copper0(1:2:256,:);%copper(256);
switch UIState
    case 0  % CPS
        n_CPS = evalin('base','n_CPS');
        Control(1).Parameters = {'Parameters',1,'startEvent',n_CPS};
        Control(2).Command = 'set&Run';
%         evalin('base',['Resource.Parameters.startEvent = ',num2str(n_CPS),';']);
%         evalin('base','Resource.DisplayWindow(1).Colormap = copper(256);');
        Control(2).Parameters = {'DisplayWindow',1,'colormap',copper0};
        Resource.DisplayWindow(1).Colormap = copper(256);
        Resource.Parameters.startEvent = n_CPS;
        
%         PData(1).Size(2) = ceil((127*Trans.spacing)/PData(1).PDelta(1));;
%         PData(1).Origin = [-Trans.spacing*127/2,0,P.startDepth];
%         PData(1).Region = computeRegions(PData(1));
%         
%         Control(3).Command = 'update&Run';
%         Control(3).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
%         
%         Control(4).Command = 'set&Run';
%         Control(5).Command = 'set&Run';
%         ScrnSize = get(0,'ScreenSize');
%         DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
%         DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
%         Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%                                               DwWidth, DwHeight];
%         Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
%         Control(4).Parameters ={'DisplayWindow',1,'Position',Resource.DisplayWindow(1).Position};
%         Control(5).Parameters ={'DisplayWindow',1,'ReferencePt',Resource.DisplayWindow(1).ReferencePt};
        Control(3).Command = 'set&Run';
        Control(3).Parameters = {'Process',1,'display',1};
        Control(4).Command = 'set&Run';
        Control(4).Parameters = {'Process',1,'mappingMode','full'};
        Process(1).display = 1;
        Process(1).mappingMethod = 'full';
        assignin('base','Process',Process);
        assignin('base','Resource',Resource)
        assignin('base','Control',Control);
%        assignin('base','PData',PData);
        assignin('base', 'action', 'displayChange');
        
%         
%         start_recording = 0;
%         assignin('base','start_recording', start_recording);
   case 1  % CPS +trigger out and cavitation detection
        n_treat = evalin('base','n_treat');
        Control(1).Parameters = {'Parameters',1,'startEvent',n_treat};
        Control(2).Command = 'set&Run';
        Resource.DisplayWindow(1).Colormap = copper(256);
%         evalin('base','Resource.DisplayWindow(1).Colormap = copper(256);');
        Control(2).Parameters = {'DisplayWindow',1,'colormap',coppercolormap};
        Resource.Parameters.startEvent = n_treat;
%         evalin('base',['Resource.Parameters.startEvent = ',num2str(n_treat),';']);
%         if kk~=1
%         kk = min(find(TX(4).apod))-1;
%         else
%             kk = min(find(TX(4).apod))
%         end

%         PData(1).Size(2) = ceil((max(find(TX(end).Apod))-min(find(TX(end).Apod)))/PData(1).PDelta(1));
%         PData(1).Origin = [Trans.ElementPos(min(find(TX(end).Apod)),1)/(Resource.Parameters.speedOfSound/Trans.frequency*1e-3),0,P.startDepth];
%         PData(1).Region = computeRegions(PData(1));
%         
%         Control(3).Command = 'update&Run';
%         Control(3).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
%         
%         Control(4).Command = 'set&Run';
%         Control(5).Command = 'set&Run';
%         ScrnSize = get(0,'ScreenSize');
%         DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
%         DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
%         Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%                                               DwWidth, DwHeight];
%         Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
%         Control(4).Parameters ={'DisplayWindow',1,'Position',Resource.DisplayWindow(1).Position};
%         Control(5).Parameters ={'DisplayWindow',1,'ReferencePt',Resource.DisplayWindow(1).ReferencePt};
        
        Control(3).Command = 'set&Run';
        Control(3).Parameters = {'Process',1,'display',0};
        Control(4).Command = 'set&Run';
        Control(4).Parameters = {'Process',1,'mappingMode','lowerHalf'};
        Process(1).display = 0;
        Process(1).mappingMethod = 'lowerHalf';
        assignin('base','Process',Process);
        assignin('base','Resource',Resource)
        assignin('base','Control',Control);
%        assignin('base','PData',PData);
        assignin('base', 'action', 'displayChange');
%         
%         start_recording = 1;
%         assignin('base','start_recording', start_recording);     
end
return
%-UI#3Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LeftEleCallback
    Trans = evalin('base','Trans');
    TX = evalin('base','TX');
    mmi = find(TX(end).Apod);
    TX(end).Apod(1:min(mmi(:))) = 1;
    ind = find(Trans.ElementPos(:,1)>UIValue);
    TX(end).Apod(1:min(ind)-1) = 0;
    assignin('base','TX',TX);
    % Set Control.Command
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);

%     bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
%     HandleImg = handle(bmodeFigHandle.Children);
%     markHandle = evalin('base','markHandle');
% if isempty(markHandle)
%     figure( HandleImg);
%     PData = evalin('base','PData');
%     P = evalin('base','P');
%     Resource = evalin('base','Resource');
%     x_plot = Trans.ElementPos(UIValue);xx = x_plot*ones(1,PData(1).Size(1));
%     yy = linspace(P.startDepth*Resource.Parameters.speedOfSound/Trans.frequency*1e-3,P.startDepth*Resource.Parameters.speedOfSound/Trans.frequency*1e-3,PData(1).Size(1));
%      hold on
%     markHandle = plot(xx,yy,'--g');
%      hold off 
% end    
%      assignin('base','markHandle',markHandle);
return
%LeftEleCallback


%RightEleCallback
    Trans = evalin('base','Trans');
    TX = evalin('base','TX');
    mma = find(TX(end).Apod);
    TX(end).Apod(max(mma(:)):end) = 1;
    ind = find(Trans.ElementPos(:,1)<UIValue);
    TX(end).Apod(max(ind)+1:end) = 0;
    assignin('base','TX',TX);
    % Set Control.Command
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);  
%     bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
%     HandleImg = handle(bmodeFigHandle.Children);
%     markHandle = evalin('base','markHandle');
% if isempty(markHandle)
%     figure(HandleImg);
%     PData = evalin('base','PData');
%     P = evalin('base','P');
%     Resource = evalin('base','Resource');
%     x_plot = Trans.ElementPos(UIValue);xx = x_plot*ones(1,PData(1).Size(1));
%     yy = linspace(P.startDepth*Resource.Parameters.speedOfSound/Trans.frequency*1e-3,P.startDepth*Resource.Parameters.speedOfSound/Trans.frequency*1e-3,PData(1).Size(1));
%      hold on
%     markHandle = plot(xx,yy,'--g');
%      hold off 
% end      
%      assignin('base','markHandle',markHandle);
return
%RightEleCallback
