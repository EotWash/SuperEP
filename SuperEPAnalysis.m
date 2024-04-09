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
files = ['Science Data\CW_2024_86a.dat'; 'Science Data\CW_2024_92b.dat'];
envFiles = ['Science Data\CW_Env_24-03-26_1118.lvm'; 'Science Data\CW_Env_24-04-01_1526.lvm'];
phse = [1 -1]; % +1 is 0 degress, -1 is 180 degrees from starting position

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
r = 3.22e-2; % Pendulum radius (m)

Msun = 1.9891e30; % Mass of sun (kg)
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Rsun = 149.6e9; % Radius from earth to sun (m)

ag = G*Msun/Rsun^2; % Acceleration towards sun (m/s^2)

CP = ((0.511/938.272)*41/93*1/41*0.05); %Cooper pair mass fraction (me/mn Ne/Nn Nv/Ne fCp);

%%

% Vector creation
inPhase = [];
outPhase = [];
unc = [];
penPhase = [];
fitTime = [];

% Long vector creation for plotting
longTorq = [];
longTim = [];
longFit = [];

% Loop through all files
for f=1:length(files(:,1))

    % Load data
    data = load(files(f,:));    

    % Parse data files into channels
    inTim = data(:,1)*3600*24; % Time vector (s)
    inDiff = data(:,2)*acCalib; % Angle vector (rad)
    
    % Calculate initial sampling frequency
    isampF = round(1/(inTim(3)-inTim(2)));
    
    % Decimate data
    inTim = decimate(inTim,isampF/decRate);
    inDiff = decimate(inDiff,isampF/decRate);
    
    % Reset sampling frequency to decimated rate
    sampF = decRate;
    
    % Trim beginning and ending of data to remove possible dead points
    startTim = 200;
    endTim = length(inTim)-100;

    tim = inTim(startTim:endTim);
    ang = inDiff(startTim:endTim);
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

    % Trim data to remove filter step response
    trim = floor(2e4*sampF);
    torqFilt = torqFilt(trim:end);
    timFilt = tim(trim:end);

    % Drift fitting to third order polynomial
    X = [0*timFilt+1 timFilt timFilt.^2 timFilt.^3];
    wt = inv(X'*X)*X'*(torqFilt);

    % Drift subtraction and vector reshaping
    torqFilt = torqFilt'-wt'*X';
    torqFilt = torqFilt';

    %% Daily Torque Fits

    % Cut size determination
    fitFreq = fDay;
    cutSize = sampF/fitFreq;
    
    % Cut based fitting
    for index = 0:floor(sampF*(timFilt(end)-timFilt(1))/cutSize)-1   

        % Cut data
        cut = torqFilt(index*cutSize+1:(index+1)*cutSize+1);    
        cutTim = timFilt(index*cutSize+1:(index+1)*cutSize+1);

        % Sync basis function and data 
        [sunMin,sunIndex]=min(abs(timSun-cutTim(1)/3600/24));
        
        % Linear least squares fitting to basis functions and offset
        X = [inSun(sunIndex:sunIndex+length(cutTim)-1) outSun(sunIndex:sunIndex+length(cutTim)-1)...
            ones(length(cut),1)];    
        w = inv(X'*X)*X'*cut;

        % Uncertainty calculation from residuals
        u = std(cut'-w'*X');

        % Vector appending
        inPhase = [inPhase w(1)];
        outPhase = [outPhase w(2)];
        unc = [unc u];
        penPhase = [penPhase phse(f)];
        fitTime = [fitTime mean(cutTim)];
        longTorq = [longTorq; cut];
        longTim = [longTim; cutTim];
        longFit = [longFit w'*X'];

    end
    % Add NaNs to plot gaps instead of lines between data sets
    longTorq = [longTorq; NaN];
    longTim = [longTim; cutTim(end)];
    longFit = [longFit NaN];

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
disp([' '])
disp(['Eta: ' num2str(etaErr*2)])
disp(['EtaCP: ' num2str(etaCPErr*2)])
%% Figures

% In-phase fits vs time
figure(1)
errorbar((fitTime-fitTime(1))/3600/24, 1e15*inPhase,1e15*unc,'.','MarkerSize',16,'LineWidth',1.5)
hold on
plot([5.5 5.5], [-250 250],'k--','LineWidth',1.5)
text(2.5,-175, '$0^{\circ}$','Interpreter', 'latex','FontSize',20)
text(8.5,-175, '$180^{\circ}$','Interpreter', 'latex','FontSize',20)
hold off
ylabel('Torque (fNm)','Interpreter', 'latex')
grid on
xlabel('Time (Days)','Interpreter', 'latex')
ylim([-250 250])
set(gca,'FontSize',16);
set(gcf, 'Position',  [50, 150, 900, 700])

% Time series of torque, fit, and in-phase function in arb. units
figure(2)
subplot(4,1,[1 3])
l=plot((longTim-longTim(1))/3600/24,1e12*(longTorq),(longTim-longTim(1))/3600/24,1e12*(longFit),timSun-timSun(1), 0.2*inSun);
legend('Observed Torque','Fit','In-Phase Function','Interpreter', 'latex','Location','southeast')
grid on
ylabel('Torque (pNm)','Interpreter', 'latex')
xlim([0 12])
ylim([-1 1])
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
set(gca,'xticklabel',[])
set(gca,'xtick',linspace(1,12,12))
subplot(4,1,4)
l=plot((longTim-longTim(1))/3600/24,1e12*(longTorq-longFit'));
legend('off')
grid on
ylabel('Residual (pNm)','Interpreter', 'latex')
xlabel('Time (Days)','Interpreter', 'latex')
xlim([0 12])
ylim([-0.4 0.4])
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
set(gca,'xtick',linspace(1,12,12))
ax = gca;
labels = string(ax.XAxis.TickLabels); 
labels(1:2:end) = nan; 
ax.XAxis.TickLabels = labels; 
set(gcf, 'Position',  [50, 100, 1500, 700])


%% Print figures

if(false)
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