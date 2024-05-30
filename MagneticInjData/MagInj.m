% warmdata = load('Warm/CW_2022_18a.dat');
% colddata = load('New/CW_2023_221a.dat');
% supdata = load('New/CW_2023_237e.dat');
% colddata = load('ColdNoGas/CW_2022_25a.dat');

warmdata = load('Warm/CW_2024_57a.dat');
colddata = load('ColdNoGas/CW_2024_58a.dat');
supdata = load('Supercon/CW_2024_150b.dat');

injData = load('Supercon/CW_Env_24-05-29_1054.lvm');
acCalib = 14e-6;

%% Room Temperature

f0 = 4.955e-3;
I = 8.574e-5;
w0=2*pi*f0;
kappa = I*w0^2;

inTim = warmdata(:,1)*3600*24;
inDiff = warmdata(:,2)*acCalib;

sampF = 1/(inTim(3)-inTim(2));
startInd = 10;
endInd = 4.4e3*sampF;
% startInd = 1e3*sampF;
% endInd = length(inTim)-10;

tim = inTim(startInd:endInd);
diff = inDiff(startInd:endInd);

tim = tim-tim(1);

Q=2e3;

filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

torqTim = kappa*lsim(filt, diff, tim);

torqTim = torqTim(500:end);

nAv = 1;
[torqW, FW] = asd2(torqTim, 1/sampF, nAv, 1, @hann);

%% Cold No Gas

f0 = 5.034e-3;
I = 8.574e-5;
w0=2*pi*f0;
kappa = I*w0^2;

inTim = colddata(:,1)*3600*24;
inDiff = colddata(:,2)*acCalib;

sampF = 1/(inTim(3)-inTim(2));
% startInd = 10;
endInd = 1.6e4*sampF;
startInd = 100;
% endInd = length(inTim);

tim = inTim(startInd:endInd);
diff = inDiff(startInd:endInd);

tim = tim-tim(1);

Q=2e3;

filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

torqTim = kappa*lsim(filt, diff, tim);

torqTim = torqTim(500:end);

nAv = 1;
[torqC, FC] = asd2(torqTim, 1/sampF, nAv, 1, @hann);


%% Cold With Gas

f0 = 5.034e-3;
I = 8.574e-5;
w0=2*pi*f0;
kappa = I*w0^2;

inTim = supdata(:,1)*3600*24;
inDiff = supdata(:,2)*acCalib;

sampF = 1/(inTim(3)-inTim(2));
startInd = 10;
% endInd = 1.6e4*sampF;
% startInd = 100;
endInd = length(inTim);

tim = inTim(startInd:endInd);
diff = inDiff(startInd:endInd);

tim = tim-tim(1);

Q=2e3;

filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

torqTim = kappa*lsim(filt, diff, tim);

torqTim = torqTim(500:end);

nAv = 11;
[torqS, FS] = asd2(torqTim, 1/sampF, nAv, 1, @hann);

%% Magnetic Coupling Measure

f0 = 5.034e-3;
I = 8.574e-5;
w0=2*pi*f0;
kappa = I*w0^2;

inMag = injData(:,9);
inTimM = injData(:,1);

inTim = supdata(:,1)*3600*24;
inDiff = supdata(:,2)*acCalib;

sampF = 1/(inTim(3)-inTim(2));
% startInd = 100;
% endInd = length(inTimM);

tim = inTim(startInd:endInd);
diff = inDiff(startInd:endInd);
mag = inMag(startInd:endInd);
timM = inTimM(startInd:endInd);

tim = tim-tim(1);

Q=2e3;

filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

torqTim = kappa*lsim(filt, diff, tim);

torqTim = torqTim(500:end);
tim = tim(500:end);
mag = mag(500:end);
timM = timM(500:end);

fitFreq = 7e-3;
cutSize = 5*sampF/fitFreq;

C = [];
S = [];

CM = [];
SM = [];

nAv = 1;
[AM, FM] = asd2(mag, (timM(2)-timM(1)), nAv, 1, @hann);

for index = 0:floor(length(tim)/cutSize)-1
   
    cut = torqTim(index*cutSize+1:(index+1)*cutSize+1);
    cutTim = tim(index*cutSize+1:(index+1)*cutSize+1);
    magCut = mag(index*cutSize+1:(index+1)*cutSize+1);    
    cutTimM = timM(index*cutSize+1:(index+1)*cutSize+1);
    
    X = [ones(length(cut),1) cutTim...
        cos(2*pi*fitFreq*cutTim) sin(2*pi*fitFreq*cutTim)];
    XM = [ones(length(cutTimM),1) cutTimM...
        cos(2*pi*fitFreq*cutTimM) sin(2*pi*fitFreq*cutTimM)];
    
    
    w = inv(X'*X)*X'*cut;
    wM = inv(XM'*XM)*XM'*magCut;
    
    C = [C w(3)];
    S = [S w(4)];
    
    CM = [CM wM(3)];
    SM = [SM wM(4)];
    
end

torqAmp = sqrt(mean(C)^2+mean(S)^2);
torqAmpErr = 1/torqAmp*sqrt(mean(C)^2*std(C)^2+mean(S)^2*std(S)^2)/sqrt(length(C));

magAmp = sqrt(mean(CM)^2+mean(SM)^2);
magAmpErr = 1/magAmp*sqrt(mean(CM)^2*std(CM)^2+mean(SM)^2*std(SM)^2)/sqrt(length(CM));

disp(['Magnetic Coupling ' num2str(torqAmp/magAmp*1e15) ' +- ' num2str(torqAmpErr/magAmp*1e15) ' fN m /Gauss']);
disp(['Approximate Systematic ' num2str(1e-7*torqAmp/magAmp*1e15) ' +- ' num2str(1e-7*torqAmpErr/magAmp*1e15) ' fN m']);
%%

fig=figure(1);
l=semilogy(FW*1e3,torqW, FC*1e3, torqC, FS*1e3, torqS,...
    [7 7], [1e-15 1e-14],'k--', [9 9], [1e-15 1e-13],'k--', [5 5], [1e-15 1e-13],'k--', [18 18], [1e-15 1e-13],'k--');
t = text(5,0.9e-15,'Resonance','Interpreter', 'latex');
set(t,'Rotation',-45);
set(t,'FontSize',16);
t = text(7,0.9e-15,'Exterior 1$\omega$','Interpreter', 'latex');
set(t,'Rotation',-45);
set(t,'FontSize',16);
t = text(9,0.9e-15,'Interior 1$\omega$','Interpreter', 'latex');
set(t,'Rotation',-45);
set(t,'FontSize',16);
t = text(18,0.9e-15,'Interior 2$\omega$','Interpreter', 'latex');
set(t,'Rotation',-45);
set(t,'FontSize',16);
xlabel('Frequency (mHz)','Interpreter', 'latex')
ylabel('Torque (N m/$\sqrt{Hz}$)','Interpreter', 'latex')
legend('Room Temperature','Cyrogenic No Exchange Gas','Cyrogenic With Exchange Gas','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
xlim([0.5 25])
ylim([1e-19 1e-6])
grid on

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'ColdWash_MagInj.pdf','-dpdf','-r1200')