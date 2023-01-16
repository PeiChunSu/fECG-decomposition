function y = step0_preprocess(x,fs)
HF_CUT = 100;% high cut frequency
[b_lp,a_lp] = butter(5,HF_CUT/(fs/2),'low');

wo = 50/(fs/2);  bw = wo/35;
[b_notch,a_notch] = iirnotch(wo,bw); %notch filter1

wo2 = 60/(fs/2);  bw2 = wo2/35;
[b_notch2,a_notch2] = iirnotch(wo2,bw2); %notch filter2

x = filtfilt(b_lp,a_lp,x);
x = filtfilt(b_notch,a_notch,x);
x = filtfilt(b_notch2,a_notch2,x);
%% remove the trend
y = ECG_detrend(x, round(fs*0.2), round(fs*0.6)); 
% preserve morphology with short window size = 200 and long window size = 600