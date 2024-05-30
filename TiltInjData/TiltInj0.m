warning('off','all')

% 0 degrees

supdata = load('CW_2024_141a.dat');
injData = load('CW_Env_24-05-20_1130.lvm');

acCalib = 14e-6;

%% Magnetic Coupling Measure

f0 = 4.94e-3;
I = 8.574e-5;
w0=2*pi*f0;
kappa = I*w0^2;

decRate = 1; % Decimate rate (Hz)

% Parse data files into channels
inTim = supdata(:,1)*3600*24; % Time vector (s)
inDiff = supdata(:,2)*acCalib; % Angle vector (rad)
inETim = injData(:,1); % Environmental time vector (s)
inTiltX = injData(:,2); % Tilt-X (rad)
inTiltY = injData(:,3); % Tilt-Y (rad) 
inICH = injData(:,4); % Inner cold head temperature (K) 
inISh = injData(:,5); % Inner shield temperature (K) 
inOCH = injData(:,6); % Outer cold head temperature (K) 
inOSh = injData(:,7); % Outer shield temperature (K) 
inFlange = 3.57*injData(:,8)+40+273.15; % Top flange temperature (K)
inWT = injData(:,9); % Water temp
inWP = injData(:,10); % Water pressure
inMagX = injData(:,11); % Magnetic field X
inMagY = injData(:,12); % Magnetic field 
inETim = inETim + inTim(1); % Add start time


% Calculate initial sampling frequency
isampF = round(1/(inTim(3)-inTim(2)));
esampF = round(1/(inETim(2)-inETim(1)));

% Decimate data
inTim = decimate(inTim,isampF/decRate);
inETim = decimate(inETim,isampF/decRate);
inDiff = decimate(inDiff,isampF/decRate);
inTiltX = decimate(inTiltX,esampF/decRate);
inTiltY = decimate(inTiltY,esampF/decRate);
inICH = decimate(inICH,esampF/decRate);
inOCH = decimate(inOCH,esampF/decRate);
inOSh = decimate(inOSh,esampF/decRate);
inISh = decimate(inISh,esampF/decRate);
inWP = decimate(inWP,esampF/decRate);
inWT = decimate(inWT,esampF/decRate);
inMagX = decimate(inMagX,esampF/decRate);
inMagY = decimate(inMagY,esampF/decRate);
inFlange = decimate(inFlange,esampF/decRate);

% Reset sampling frequency to decimated rate
sampF = decRate;

%% Y-tilt

startInd = 3e3;
endInd = 1.4e4;

tim = inTim(startInd:endInd);
ang = inDiff(startInd:endInd);
tiltX = inTiltX(startInd:endInd);
tiltY = inTiltY(startInd:endInd);

tim = tim-tim(1);

trim = 10;
tim = tim(trim:end);
tiltX = tiltX(trim:end);
tiltY = tiltY(trim:end);
ang = ang(trim:end);

angPlot = movmean(ang-mean(ang),201);

angPlot = angPlot(200:end-200);
tiltYPlot = tiltY(200:end-200);
X = [ones(length(tiltYPlot),1) tiltYPlot];
y = angPlot;
w = inv(X'*X)*X'*y;
coupy = w(2);
disp(['Tilt-Y Coupling ' num2str(coupy) ' rad/rad']);

fig=figure(1);
l=plot(tim,tiltX,tim,tiltY,tim,movmean(ang-mean(ang),201));
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Angle (rad)','Interpreter', 'latex')
legend('Tilt-X','Tilt-Y','Angle','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on

fig=figure(2);
l=plot(tiltYPlot,angPlot,'.',tiltYPlot,w'*X');
xlabel('Tilt-Y (rad)','Interpreter', 'latex')
ylabel('Angle (rad)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on

%% X- Tilt

startInd = 1.6e5;
endInd = length(inTim)-100;
% endInd = 2e4;

tim = inTim(startInd:endInd);
ang = inDiff(startInd:endInd);
tiltX = inTiltX(startInd:endInd);
tiltY = inTiltY(startInd:endInd);

tim = tim-tim(1);

trim = 10;
tim = tim(trim:end);
tiltX = tiltX(trim:end);
tiltY = tiltY(trim:end);
ang = ang(trim:end);

%%
angPlot = movmean(ang-mean(ang),201);

angPlot = angPlot(200:end-200);
tiltXPlot = tiltX(200:end-200);
tiltYPlot = tiltY(200:end-200);
angPlot = angPlot-coupy*tiltYPlot;

X = [ones(length(tiltXPlot),1) tiltXPlot];
y = angPlot;
w = inv(X'*X)*X'*y;
coupx = w(2);

disp(['Tilt-X Coupling ' num2str(coupx) ' rad/rad']);

fig=figure(3);
l=plot(tim,tiltX,tim,tiltY,tim,movmean(ang-mean(ang),201));
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Angle (rad)','Interpreter', 'latex')
legend('Tilt-X','Tilt-Y','Angle','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on

fig=figure(4);
l=plot(tiltXPlot,angPlot,'.',tiltXPlot,w'*X');
xlabel('Tilt-X (rad)','Interpreter', 'latex')
ylabel('Angle (rad)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
grid on

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'ColdWash_MagInj.pdf','-dpdf','-r1200')