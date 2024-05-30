warning('off','all')

% Artifical signal injection control
inj = false;
injAmp = 100e-15;

% Loading in basis funtions outputted from sunVect.py
rawSun=load('sunVectMin.out');
timSun=rawSun(:,1);
inSun=detrend(rawSun(:,2),'linear');
outSun=detrend(rawSun(:,3),'linear');

% File paths and pendulum orientations
files = ["Science Data\CW_2024_86a.dat" "Science Data\CW_2024_92b.dat" "Science Data\CW_2024_99a.dat" "Science Data\CW_2024_108a.dat"];
envFiles = ["Science Data\CW_Env_24-03-26_1118.lvm" "Science Data\CW_Env_24-04-01_1526.lvm" ...
    "Science Data\CW_Env_24-04-08_1312.lvm" "Science Data\CW_Env_24-04-17_1135.lvm"];
phse = [1 -1 1 -1]; % +1 is 0 degress, -1 is 180 degrees from starting position

%%


acCalib = 14e-6; % Autocollimator calibration (rad/pixel)

f0 = 5.07e-3; % Pendulum resonance (Hz)
I = 8.574e-5; % Pendulum moment of inertia (kg m^2)
w0 = 2*pi*f0; % Resonance in angular frequency units (rad/s)
kappa = I*w0^2; % Spring constant (N m/rad)
Q = 16198; % Quality factor
fDay = 1/(3600*24); % Once per day frequency (Hz)

decRate = 1/60; % Decimate rate (Hz)

m = 37.92e-3; % Pendulum masses (kg)
r = 2.27e-2; % Pendulum radius (m)

Msun = 1.9891e30; % Mass of sun (kg)
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Rsun = 149.6e9; % Radius from earth to sun (m)

ag = G*Msun/Rsun^2; % Acceleration towards sun (m/s^2)

CP = ((0.511/938.272)*41/93*1/41*0.05); % Cooper pair mass fraction (me/mn Ne/Nn Nv/Ne fCp);

tempCoup = 4.6086e-11; % Inner cold-head temperature coupling (N m/K)
txCoup = -0.62704*kappa; % Tilt-X coupling (N m/rad)
tyCoup = 0.40461*kappa; % Tilt-Y coupling (N m/rad)

%%

% Vector creation
inPhase = [];
outPhase = [];
unc = [];
penPhase = [];
fitTime = [];

% Long vector creation for plotting
longTorq = [];
longICH = [];
longTim = [];
longFit = [];

% Loop through all files
for f=1:length(files)

    % Load data
    data = load(files(f));
    envData = load(envFiles(f)); 

    % Parse data files into channels
    inTim = data(:,1)*3600*24; % Time vector (s)
    inDiff = data(:,2)*acCalib; % Angle vector (rad)
    inETim = envData(:,1); % Environmental time vector (s)
    inTiltX = envData(:,2); % Tilt-X (rad)
    inTiltY = envData(:,3); % Tilt-Y (rad) 
    inICH = envData(:,4); % Inner cold head temperature (K) 
    inISh = envData(:,5); % Inner shield temperature (K) 
    inOCH = envData(:,6); % Outer cold head temperature (K) 
    inOSh = envData(:,7); % Outer shield temperature (K) 
    inFlange = 3.57*envData(:,8)+40+273.15; % Top flange temperature (K)
    inWT = envData(:,9); % Water temp
    inWP = envData(:,10); % Water pressure
    inMagX = envData(:,11); % Magnetic field X
    inMagY = envData(:,12); % Magnetic field 
    inETim = inETim + inTim(1); % Add start time
    
    
    % Calculate initial sampling frequency
    isampF = round(1/(inTim(3)-inTim(2)));
    esampF = round(1/(inETim(2)-inETim(1)));

    % Decimate data
    inTim = decimate(inTim,isampF/decRate);
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
    
    % Trim beginning and ending of data to remove possible dead points
    startTim = 200;
    endTim = length(inTim)-100;

    tim = inTim(startTim:endTim);
    ang = inDiff(startTim:endTim);
    tiltX = inTiltX(startTim:endTim);
    tiltY = inTiltY(startTim:endTim);
    ICH = inICH(startTim:endTim);
    ISh = inISh(startTim:endTim);
    OCH = inOCH(startTim:endTim);
    OSh = inOSh(startTim:endTim);
    WP = inWP(startTim:endTim);
    WT = inWT(startTim:endTim);
    magX = inMagX(startTim:endTim);
    magY = inMagY(startTim:endTim);
    flange = inFlange(startTim:endTim);
    %% Torque Calculations

    % Calculate torque
    torq = kappa*ang; % Torque vector (N m)
    
    % Signal injection
    if (inj)
        % Sync basis function and data
        [sunMin,sunIndex]=min(abs(timSun-tim(1)/3600/24));
        torq = torq + phse(f)*injAmp*inSun(sunIndex:sunIndex+length(tim)-1);
    end

    % Low pass filter to remove resonance and high frequency noise
    [b,a] = butter(2,2*1e-4/sampF,'low');
    torqFilt = filter(b,a,torq);
    tXFilt = filter(b,a,tiltX);
    tYFilt = filter(b,a,tiltY);
    ICHFilt = filter(b,a,ICH);
    IShFilt = filter(b,a,ISh);
    OCHFilt = filter(b,a,OCH);
    OShFilt = filter(b,a,OSh);
    WPFilt = filter(b,a,WP);
    WTFilt = filter(b,a,WT);
    magXFilt = filter(b,a,magX);
    magYFilt = filter(b,a,magY);
    flFilt = filter(b,a,flange);

    % Trim data to remove filter step response
    trim = floor(2e4*sampF);
    torqFilt = torqFilt(trim:end);
    timFilt = tim(trim:end); 
    tXFilt = tXFilt(trim:end);
    tYFilt = tYFilt(trim:end);
    ICHFilt =ICHFilt(trim:end);
    IShFilt = IShFilt(trim:end);
    OCHFilt = OCHFilt(trim:end);
    OShFilt = OShFilt(trim:end);
    WPFilt = WPFilt(trim:end);
    WTFilt = WTFilt(trim:end);
    magXFilt = magXFilt(trim:end);
    magYFilt = magYFilt(trim:end);
    flFilt = flFilt(trim:end);

    % Drift fitting to third order polynomial
    X = [0*timFilt+1 timFilt timFilt.^2 timFilt.^3];
    wt = inv(X'*X)*X'*(torqFilt);
    w = inv(X'*X)*X'*(ICHFilt);
    wx = inv(X'*X)*X'*(tXFilt);
    wy = inv(X'*X)*X'*(tYFilt);

    % Drift subtraction and vector reshaping
    torqFilt = torqFilt'-wt'*X';
    torqFilt = torqFilt'-mean(torqFilt);
    ICHFilt = ICHFilt'-w'*X';
    ICHFilt = ICHFilt';
    tXFilt = tXFilt'-wx'*X';
    tXFilt = tXFilt';
    tYFilt = tYFilt'-wy'*X';
    tYFilt = tYFilt';


    %% Daily Torque Fits

    % Cut size determination
    fitFreq = fDay;
    cutSize = sampF/fitFreq;
    
    % Cut based fitting
    for index = 0:floor(sampF*(timFilt(end)-timFilt(1))/cutSize)-1   

        % Cut data
        cut = torqFilt(index*cutSize+1:(index+1)*cutSize+1);    
        cutTim = timFilt(index*cutSize+1:(index+1)*cutSize+1);
        cutICH = ICHFilt(index*cutSize+1:(index+1)*cutSize+1);
        cutTX = tXFilt(index*cutSize+1:(index+1)*cutSize+1);
        cutTY = tYFilt(index*cutSize+1:(index+1)*cutSize+1);
        
        % Temperature correction
        cutCor = cut - tempCoup*cutICH -txCoup*cutTX -tyCoup*cutTY;
%         cutCor = cut;

        % Sync basis function and data 
        [sunMin,sunIndex]=min(abs(timSun-cutTim(1)/3600/24));
        
        % Linear least squares fitting to basis functions and offset
        X = [inSun(sunIndex:sunIndex+length(cutTim)-1) outSun(sunIndex:sunIndex+length(cutTim)-1)...
            ones(length(cut),1)];    
        w = inv(X'*X)*X'*cutCor;

        % Uncertainty calculation from residuals
        u = std(cutCor'-w'*X');

        % Vector appending
        inPhase = [inPhase w(1)];
        outPhase = [outPhase w(2)];
        unc = [unc u];
        penPhase = [penPhase phse(f)];
        fitTime = [fitTime mean(cutTim)];
        longTorq = [longTorq; cutCor];
        longICH = [longICH; cutICH];
        longTim = [longTim; cutTim];
        longFit = [longFit w'*X'];
        
        % Add NaNs to plot gaps instead of lines between data sets
        longTorq = [longTorq; NaN];
        longICH = [longICH; NaN];
        longTim = [longTim; cutTim(end)];
        longFit = [longFit NaN];

    end


end
%%

% Uncertainty weighted average of in-phase amplitude
mu = sum(penPhase.*inPhase./(unc.^2))/sum(1./(unc.^2));
muErr = sqrt(1/sum(1./(unc.^2)));
% Uncertainty weighted average of out-of-phase amplitude
muOut = sum(penPhase.*outPhase./(unc.^2))/sum(1./(unc.^2));

% Conversion to Eotvos parameter for tungsten
eta = mu/(r*m*ag);
etaErr = muErr/(r*m*ag);

% Conversion to Eotvos parameter for Cooper pairs
etaCP= eta/CP;
etaCPErr= etaErr/CP;

% Display results
disp(['Torque: ' num2str(mu*1e15) ' +- ' num2str(muErr*1e15)])
disp(['Torque Out: ' num2str(muOut*1e15) ' +- ' num2str(muErr*1e15)])
disp([' '])
disp(['Eta: ' num2str(etaErr*2)])
disp(['EtaCP: ' num2str(etaCPErr*2)])
%% Figures

Nday = floor(1.05*(longTim(end)-longTim(1))/24/3600);

% In-phase fits vs time
figure(1)
errorbar((fitTime-fitTime(1))/3600/24, 1e15*inPhase,1e15*unc,'.','MarkerSize',16,'LineWidth',1.5)
hold on
plot([5.5 5.5], [-250 250],'k--','LineWidth',1.5)
plot([12 12], [-250 250],'k--','LineWidth',1.5)
plot([21 21], [-250 250],'k--','LineWidth',1.5)
text(2.5,-175, '$0^{\circ}$','Interpreter', 'latex','FontSize',20)
text(8.5,-175, '$180^{\circ}$','Interpreter', 'latex','FontSize',20)
text(16,-175, '$0^{\circ}$','Interpreter', 'latex','FontSize',20)
text(23,-175, '$180^{\circ}$','Interpreter', 'latex','FontSize',20)
hold off
ylabel('Torque (fNm)','Interpreter', 'latex')
grid on
xlabel('Time (Days)','Interpreter', 'latex')
ylim([-250 250])
xlim([-1 Nday])
set(gca,'FontSize',16);
set(gcf, 'Position',  [50, 150, 900, 700])

% Time series of torque, fit, and in-phase function in arb. units
figure(2)
subplot(4,1,[1 3])
l=plot((longTim-longTim(1))/3600/24,1e12*(longTorq),(longTim-longTim(1))/3600/24,1e12*(longFit),timSun-timSun(1), 0.2*inSun);
legend('Observed Torque','Fit','In-Phase Function','Cold Head Temp', 'Raw Torque','Interpreter', 'latex','Location','southeast')
grid on
ylabel('Torque (pNm)','Interpreter', 'latex')
xlim([-1 Nday])
ylim([-0.5 0.5])
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
set(gca,'xticklabel',[])
set(gca,'xtick',linspace(0,Nday,Nday+1))
set(gca,'ytick',linspace(-0.5,0.5,9))
subplot(4,1,4)
l=plot((longTim-longTim(1))/3600/24,1e12*(longTorq-longFit'));
legend('off')
grid on
ylabel('Residual (pNm)','Interpreter', 'latex')
xlabel('Time (Days)','Interpreter', 'latex')
xlim([-1 Nday])
ylim([-0.41 0.41])
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
set(gca,'xtick',linspace(0,Nday,Nday+1))
ax = gca;
labels = string(ax.XAxis.TickLabels); 
labels(2:2:end) = nan; 
ax.XAxis.TickLabels = labels; 
set(gcf, 'Position',  [50, 100, 1500, 700])


%% Print figures

if(true)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'SuperEP_Fit.pdf','-dpdf','-r1200')
    
    fig2=figure(2);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'SuperEP_TimeSeries.pdf','-dpdf','-r1200')
end