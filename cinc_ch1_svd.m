
clear; close all ;

DBfolder ='D:\duke box\fECG_project\Code Archive3\cinc2013\set-a-text';
addpath(genpath('D:\duke box\fECG_project\Code Archive3\Code Archive3'))


total_F1 = zeros(75,4);
total_MAE = zeros(75,4);
total_PPV = zeros(75,4);
total_SE = zeros(75,4);
total_time = zeros(75,4);

final_mbeats = {};
final_fbeats = {};
final_aECG = {};
final_Om = {};
final_I2_orig = {};
final_Of = {};
alpha_c = [];


%% ==================== SQI paremeters set up================================
opt = struct(...
    'SIZE_WIND',4,... % define the window for the bSQI check on the ECG
    'LG_MED',0,... % take the median SQI using X nearby values,  so if LG_MED = 3, we take the median of the 3 prior and 3 posterior windows
    'REG_WIN',1,... % how frequently to check the SQI for switching - i.e., if REG_WIN = 1, then we check the signals every second to switch
    'THR',0.150,... % the width, in seconds, used when comparing peaks in the F1 based ECG SQI
    'SQI_THR',0.8,... % the SQI threshold - we switch signals if SQI < this value
    'USE_PACING',1,... % flag turning on/off the pacing detection/correction
    'ABPMethod','wabp',... % ABP peak detection method (wabp, delineator)
    'SIMPLEMODE',0,... % simple mode only uses the first ABP and ECG signal, and ignores all others
    'DELAYALG', 'map',... % algorithm used to determine the delay between the ABP and the ECG signal
    ... % jqrs parameters - the custom peak detector implemented herein
    'JQRS_THRESH', 0.3,... % energy threshold for defining peaks
    'JQRS_REFRAC', 0.25,... % refractory period in seconds
    'JQRS_INTWIN_SZ', 7,...
    'JQRS_WINDOW', 15);


%% ============================== deshape parameters setup ==================

basicTF.win = 500;
basicTF.hop = 20;
basicTF.fs = 100;
basicTF.fr = 0.02;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Flattop'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
%advTF.win_type = 'Hamming'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-6;
advTF.HighFreq = 10/100;
advTF.LowFreq = 0.5/100;
advTF.lpc = 0;
cepR.g = 0.3; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;
num_tap = 1;
lam_curve = 20 ;% lambda for curve extraction (suggest: 10)
lam_beat = 50 ;% lambda for beat tracking (suggest: 20)
alphaN = 7; % # direction of linear combinations = 2*alphaN+1

%% ============================== morphology parameters setup ===============

num_med = 101;% length of median filter (suggest: 101)
num_med_l = 601;% length of long window median filter for morphology preservation (suggest: 601)
num_med_s = 201;% length of short window median filter for morphology preservation (suggest: 201)
num_nonlocal = 20;% number of waveform for nonlocal median
OptimalShrinkageOpt = 'op2'; % Optimal shrinkage norm options: 'fro', 'nuc', 'op'


%% ============================== preprocess parameters setup ===============
fs = 1000; %#sampling rate
HF_CUT = 100;% high cut frequency
[b_lp,a_lp] = butter(5,HF_CUT/(fs/2),'low');
wo = 50/(1000/2);  bw = wo/35;
[b_notch,a_notch] = iirnotch(wo,bw); %notch filter
wo2 = 60/(1000/2);  bw2 = wo2/35;
[b_notch2,a_notch2] = iirnotch(wo2,bw2); %notch filter

%{
%% ============================ Loading cinc2013 Database ===================
for sub_i = 1:75
    dname = [DBfolder '/a' sprintf('%02d',sub_i) '.csv'];
    recorddata = readtable([DBfolder '/a' sprintf('%02d',sub_i) '.csv']);
    recorddata = recorddata{2:end,:};
    index = find(strcmp(recorddata, '-'));
    recorddata(index) = {'-32768'};
    recorddata = sprintf('%s ', recorddata{:});
    recorddara = sscanf(recorddata, '%f');
    recorddata = reshape(recorddara,[60000,5]);
    recorddata = recorddata(:,2:end);
    recorddata = recorddata';
    
    for ch_i = 1:4
        sig_cha = recorddata(ch_i,:);
        tmp = find(sig_cha~=-32768) ;
        sig_cha = interp1(tmp, sig_cha(tmp), [1:length(sig_cha)],'pchip','extrap') ;
        if length(tmp)<length(sig_cha) ; disp('ERROR in the data. Corrected by spline') ; end
        sig_c(sub_i,:,ch_i) = sig_cha;
    end
end
save('sig_c', 'sig_c');
load('sig_c.mat');

%% ============================ Start Counting Time =========================

disp([num_med num_nonlocal lam_curve lam_beat]) ;
TTT  = tic ;


%% =========================== Preprocess ===================================

x_up_c = {}; % collection of signals without morphology
x_up_c_real = {}; % collection of signals with morphology

cc = 0;
for (sub_i = 1:75)
    
    disp(sub_i)
    for ch_i = 1:4
        cc = cc+1;
        sig_cha_1 = sig_c(sub_i,:,ch_i);
        sig_chq_1 = filtfilt(b_lp,a_lp,sig_cha_1);
        sig_cha_1 = filtfilt(b_notch,a_notch,sig_cha_1);
        sig_cha_1 = filtfilt(b_notch2,a_notch2,sig_cha_1);
        %% remove the trend
        x0 = ECG_detrend(sig_cha_1, num_med, num_med_l, 0); % no morphology
        x0_real = ECG_detrend(sig_cha_1, num_med_s, num_med_l, 1); % preserve morphology
        %% Compute mbeats
        
        x_up_c{cc} = x0;
        x_up_c_real{cc} = x0_real;
        
    end
end

save('x_up_c_real', 'x_up_c_real');
save('x_up_c', 'x_up_c');
%}
load('x_up_c.mat');
load('x_up_c_real.mat');


for sub_i = 1:75
    sub_i
    PAIRidx = 0 ;
    FILENAME = ['/a' sprintf('%02d',sub_i)] ;
    [ann0] = load([DBfolder '/' FILENAME '.fqrs.txt']);
    
    for ch_i = 2;%1:4
        tic;
        gtbeats = ann0;
        % Read the data
        PAIRidx = PAIRidx + 1 ;
        
        x0 = x_up_c{4*(sub_i-1)+ch_i};
        x0_real = x_up_c_real{4*(sub_i-1)+ch_i};
        [x0, x0_real,mbeats, ~, ~, ~, ~, HR_ma, tfrtic, t] = mbeats_extract_real(x0,x0_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, 1);
        [Om,Om_real] = ECG_shrinkage0( x0,x0_real, mbeats, OptimalShrinkageOpt,1.5,0);
        %[Om,Om_real] = ECG_shrinkage0_color( x0,x0_real, mbeats, OptimalShrinkageOpt);   
        
        %% Remove maternal ECG waveforms
        %% polarization is applied
        %I2_orig = sqrt(1-alpha^2).*(sig_cha_1-Om1) + alpha .* (sig_cha_2-Om2);
        I2_orig_real = x0_real-Om_real;
        I2_orig = x0-Om;
        %find fECG peak locations with fusion mECG
        [~,~, fbeats, ~, ~, ~, tfrrF, HR_fe, tfrtic, t] = fbeats_extract_real(I2_orig,I2_orig_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, 0, HR_ma);
        if length(fbeats) == 0
            Om = compute_morph_nonlocal_median(x0, mbeats,num_nonlocal, 1000);
            Om_real = compute_morph_nonlocal_median(x0_real, mbeats,num_nonlocal, 1000);
            I2_orig = x0-Om; I2_orig_real = x0_real-Om_real;
            [~,~, fbeats, ~, ~, ~, tfrrF, HR_fe, tfrtic, t] = fbeats_extract_real(I2_orig,I2_orig_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, 0, HR_ma);
        end
        
        [Of,Of_real] = ECG_shrinkage0( I2_orig,I2_orig_real, fbeats, OptimalShrinkageOpt,1.5,0);
        %[Of,Of_real] = ECG_shrinkage0_color( I2_orig,I2_orig_real, fbeats, OptimalShrinkageOpt);    
        %Of_real = compute_morph_nonlocal_median(I2_orig_real, fbeats,num_nonlocal, 1000);
        
        Om_raw = x0-Of; Om_raw_real = x0_real - Of_real;
        [~, ~,mbeats, ~, ~, ~, tfrrM, HR_ma, tfrtic, t] = mbeats_extract_real(Om_raw,Om_raw_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, 1);
        [Om,Om_real] = ECG_shrinkage0( Om_raw,Om_raw_real, mbeats, OptimalShrinkageOpt,1.5,0);
        %[Om,Om_real] = ECG_shrinkage0_color( Om_raw,Om_raw_real, mbeats, OptimalShrinkageOpt);     
        I2_orig = x0-Om; I2_orig_real = x0_real-Om_real;
        
        [~,~, fbeats, ~, ~, ~, tfrrF, HR_fe, tfrtic, t] = fbeats_extract_real(I2_orig,I2_orig_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, 0, HR_ma);
        %[Of,Of_real] = ECG_shrinkage_median( I2_orig,I2_orig_real, fbeats, OptimalShrinkageOpt,num_nonlocal,0);
        
        
        
        
        
        clear Ot tempb
        %%
        
        
        mbeats = mbeats(mbeats>2000); mbeats = mbeats(mbeats<58000);
        fbeats = fbeats(fbeats>2000); fbeats = fbeats(fbeats<58000);
        gtbeats = gtbeats(gtbeats>2000); gtbeats = gtbeats(gtbeats<58000);
        
        [F1_m,MAE,PPV,SE,TP,FN,FP] = Bxb_compare(gtbeats, mbeats, 50);
        [F1_f,MAE,PPV,SE,TP,FN,FP] = Bxb_compare(gtbeats, fbeats, 50);
        
        if median(diff(fbeats))>median(diff(mbeats))
            [fbeats, mbeats] = deal(mbeats, fbeats);
            [Of_real, Om_real] = deal(Om_real, Of_real);
        end
        
        
        [F1,MAE,PPV,SE,TP,FN,FP] = Bxb_compare(gtbeats, fbeats, 50);
        
        
        F1
        total_F1(sub_i, ch_i) = F1 ;
        total_MAE(sub_i, ch_i) = MAE;
        total_PPV(sub_i, ch_i) = PPV;
        total_SE(sub_i, ch_i) = SE;
        total_TP(sub_i, ch_i) = TP;
        total_FN(sub_i, ch_i) = FN;
        total_FP(sub_i, ch_i) = FP;
        total_time(sub_i, ch_i) = toc;
        final_mbeats{sub_i,PAIRidx} = mbeats;
        final_fbeats{sub_i,PAIRidx} = fbeats;
        final_aECG{sub_i,PAIRidx} = x0_real;
        final_Om{sub_i,PAIRidx} = Om_real;
        final_Of{sub_i,PAIRidx} = Of_real;
        final_I2_orig{sub_i,PAIRidx}  = I2_orig_real;
        
        
        if F1 < 0.9
            
            %disp(F1)
            if 0
                subplot(length(alpha_ind), 1, 1) ; title(['F1 = ',num2str(F1)]) ;
                export_fig([num2str(sub_i),'_',num2str(PAIRidx)],'-transparent','-m2') ;
            end
        end
    end
end
save('cinc_shrinkage_ch1_F1_0619_color','total_F1', 'total_MAE', 'total_PPV', 'total_SE','total_TP','total_FN','total_FP','total_time','alpha_c');
%save('cinc_shrinkage_ch1_signals_0619','final_mbeats','final_fbeats','final_aECG','final_Om','final_Of','final_I2_orig');

%disp([num2str(error_count) ' errors occuredd'])



BADAndreotti2014 = [54] ;
GOODA2014 = setdiff([1:75],BADAndreotti2014) ;
BADBehar2014 = [27 32 43 50 54 59 60 63 64 68 75] ;
GOODB2014 = setdiff([1:75], BADBehar2014) ;

CH = []; for ii=1:75; [a,b]=sort(total_F1(ii,:),'descend'); CH=[CH b(1)]; end
TTP = 0 ; TFP = 0 ; TFN = 0 ;
for ii=GOODA2014; TTP=TTP+total_TP(ii,CH(ii)); TFN=TFN+total_FN(ii,CH(ii)); TFP=TFP+total_FP(ii,CH(ii)); end

QF1=[]; for ii=GOODA2014; QF1=[QF1 total_F1(ii,CH(ii))]; end
QMAE=[]; for ii=GOODA2014; QMAE=[QMAE total_MAE(ii,CH(ii))]; end
QPPV=[]; for ii=GOODA2014; QPPV=[QPPV total_PPV(ii,CH(ii))]; end
QSE=[]; for ii=GOODA2014; QSE=[QSE total_SE(ii,CH(ii))]; end

disp(['(best, No A) Total F1 all =',num2str(2*TTP./(2*TTP+TFP+TFN))]) ;
disp(['(best, No A) Mean (std) F1 all = ',num2str(mean(QF1)),' ',num2str(std(QF1))]) ;
disp(['(best, No A) Mean (std) PPV all = ',num2str(mean(QPPV)),' ',num2str(std(QPPV))]) ;
disp(['(best, No A) Mean (std) SE all = ',num2str(mean(QSE)),' ',num2str(std(QSE))]) ;
disp(['(best, No A) Mean (std) MAE all = ',num2str(mean(QMAE)),' ',num2str(std(QMAE))]) ;

disp(['(best, No A) Median (IQR) F1 all = ',num2str(median(QF1)),' ',num2str(iqr(QF1))]) ;
disp(['(best, No A) Median (IQR) PPV all = ',num2str(median(QPPV)),' ',num2str(iqr(QPPV))]) ;
disp(['(best, No A) Median (IQR) SE all = ',num2str(median(QSE)),' ',num2str(iqr(QSE))]) ;
disp(['(best, No A) Median (IQR) MAE all = ',num2str(median(QMAE)),' ',num2str(iqr(QMAE))]) ;


TTP = 0 ; TFP = 0 ; TFN = 0 ;
for ii=GOODB2014; TTP=TTP+total_TP(ii,CH(ii)); TFN=TFN+total_FN(ii,CH(ii)); TFP=TFP+total_FP(ii,CH(ii)); end

QF1=[]; for ii=GOODB2014; QF1=[QF1 total_F1(ii,CH(ii))]; end
QMAE=[]; for ii=GOODB2014; QMAE=[QMAE total_MAE(ii,CH(ii))]; end
QPPV=[]; for ii=GOODB2014; QPPV=[QPPV total_PPV(ii,CH(ii))]; end
QSE=[]; for ii=GOODB2014; QSE=[QSE total_SE(ii,CH(ii))]; end

disp(['(best, No B) Total F1 all =',num2str(2*TTP./(2*TTP+TFP+TFN))]) ;
disp(['(best, No B) Mean (std) F1 all = ',num2str(mean(QF1)),' ',num2str(std(QF1))]) ;
disp(['(best, No B) Mean (std) PPV all = ',num2str(mean(QPPV)),' ',num2str(std(QPPV))]) ;
disp(['(best, No B) Mean (std) SE all = ',num2str(mean(QSE)),' ',num2str(std(QSE))]) ;
disp(['(best, No B) Mean (std) MAE all = ',num2str(mean(QMAE)),' ',num2str(std(QMAE))]) ;

disp(['(best, No B) Median (IQR) F1 all = ',num2str(median(QF1)),' ',num2str(iqr(QF1))]) ;
disp(['(best, No B) Median (IQR) PPV all = ',num2str(median(QPPV)),' ',num2str(iqr(QPPV))]) ;
disp(['(best, No B) Median (IQR) SE all = ',num2str(median(QSE)),' ',num2str(iqr(QSE))]) ;
disp(['(best, No B) Median (IQR) MAE all = ',num2str(median(QMAE)),' ',num2str(iqr(QMAE))]) ;



%===============
fprintf('=======================\n')


fprintf(['(No A) 1 : median F1 = ',num2str(median(total_F1(GOODA2014,1))),' (IQR = ',num2str(iqr(total_F1(GOODA2014,1))),')\n']) ;
fprintf(['(No A) 2 : median F1 = ',num2str(median(total_F1(GOODA2014,2))),' (IQR = ',num2str(iqr(total_F1(GOODA2014,2))),')\n']) ;
fprintf(['(No A) 3 : median F1 = ',num2str(median(total_F1(GOODA2014,3))),' (IQR = ',num2str(iqr(total_F1(GOODA2014,3))),')\n']) ;
fprintf(['(No A) 4 : median F1 = ',num2str(median(total_F1(GOODA2014,4))),' (IQR = ',num2str(iqr(total_F1(GOODA2014,4))),')\n']) ;

fprintf(['(No A) 1 : mean F1 = ',num2str(mean(total_F1(GOODA2014,1))),' (std = ',num2str(std(total_F1(GOODA2014,1))),')\n']) ;
fprintf(['(No A) 2 : mean F1 = ',num2str(mean(total_F1(GOODA2014,2))),' (std = ',num2str(std(total_F1(GOODA2014,2))),')\n']) ;
fprintf(['(No A) 3 : mean F1 = ',num2str(mean(total_F1(GOODA2014,3))),' (std = ',num2str(std(total_F1(GOODA2014,3))),')\n']) ;
fprintf(['(No A) 4 : mean F1 = ',num2str(mean(total_F1(GOODA2014,4))),' (std = ',num2str(std(total_F1(GOODA2014,4))),')\n']) ;

fprintf(['(No B) 1 : median F1 = ',num2str(median(total_F1(GOODB2014,1))),' (IQR = ',num2str(iqr(total_F1(GOODB2014,1))),')\n']) ;
fprintf(['(No B) 2 : median F1 = ',num2str(median(total_F1(GOODB2014,2))),' (IQR = ',num2str(iqr(total_F1(GOODB2014,2))),')\n']) ;
fprintf(['(No B) 3 : median F1 = ',num2str(median(total_F1(GOODB2014,3))),' (IQR = ',num2str(iqr(total_F1(GOODB2014,3))),')\n']) ;
fprintf(['(No B) 4 : median F1 = ',num2str(median(total_F1(GOODB2014,4))),' (IQR = ',num2str(iqr(total_F1(GOODB2014,4))),')\n']) ;

fprintf(['(No B) 1 : mean F1 = ',num2str(mean(total_F1(GOODB2014,1))),' (std = ',num2str(std(total_F1(GOODB2014,1))),')\n']) ;
fprintf(['(No B) 2 : mean F1 = ',num2str(mean(total_F1(GOODB2014,2))),' (std = ',num2str(std(total_F1(GOODB2014,2))),')\n']) ;
fprintf(['(No B) 3 : mean F1 = ',num2str(mean(total_F1(GOODB2014,3))),' (std = ',num2str(std(total_F1(GOODB2014,3))),')\n']) ;
fprintf(['(No B) 4 : mean F1 = ',num2str(mean(total_F1(GOODB2014,4))),' (std = ',num2str(std(total_F1(GOODB2014,4))),')\n']) ;

fprintf(['(No A) 1 : median MAE = ',num2str(median(total_MAE(GOODA2014,1))),' (IQR = ',num2str(iqr(total_MAE(GOODA2014,1))),')\n']) ;
fprintf(['(No A) 2 : median MAE = ',num2str(median(total_MAE(GOODA2014,2))),' (IQR = ',num2str(iqr(total_MAE(GOODA2014,2))),')\n']) ;
fprintf(['(No A) 3 : median MAE = ',num2str(median(total_MAE(GOODA2014,3))),' (IQR = ',num2str(iqr(total_MAE(GOODA2014,3))),')\n']) ;
fprintf(['(No A) 4 : median MAE = ',num2str(median(total_MAE(GOODA2014,4))),' (IQR = ',num2str(iqr(total_MAE(GOODA2014,4))),')\n']) ;

fprintf(['(No A) 1 : mean MAE = ',num2str(mean(total_MAE(GOODA2014,1))),' (std = ',num2str(std(total_MAE(GOODA2014,1))),')\n']) ;
fprintf(['(No A) 2 : mean MAE = ',num2str(mean(total_MAE(GOODA2014,2))),' (std = ',num2str(std(total_MAE(GOODA2014,2))),')\n']) ;
fprintf(['(No A) 3 : mean MAE = ',num2str(mean(total_MAE(GOODA2014,3))),' (std = ',num2str(std(total_MAE(GOODA2014,3))),')\n']) ;
fprintf(['(No A) 4 : mean MAE = ',num2str(mean(total_MAE(GOODA2014,4))),' (std = ',num2str(std(total_MAE(GOODA2014,4))),')\n']) ;

F105_mean = mean(median(total_F1'))
F105_std = std(median(total_F1'))
MAE05_mean = mean(median(total_MAE'))
MAE05_std = std(median(total_MAE'))

