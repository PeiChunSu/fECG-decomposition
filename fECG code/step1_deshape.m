function tfrr = step1_deshape(x,fs,f_r,gamma,alpha,win_len,win_type)

basicTF.hop = 20;
basicTF.fs = 100;
basicTF.win = win_len*100;
basicTF.fr = f_r;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = win_type; %{'Flattop','Gaussian','Thomson','multipeak','SWCE'}; %}
%advTF.win_type = 'Hamming'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-6;
advTF.HighFreq = 10/100;
advTF.LowFreq = 0.5/100;
advTF.lpc = 0;
cepR.g = gamma; % g(cc);%0.06;
cepR.Tc=0;
cepR.alpha = alpha;
P.num_s = 1;
P.num_c = 1;


x_len = floor(length(x)/200)*200;
x = x(1:x_len);
x1 = resample(x,basicTF.fs,fs);
gg = x1; gg = gg-mean(gg);
[~, ~, ~, tfrr, ~, ~, tfrtic, t] = CFPH(gg, basicTF, advTF, cepR, P);