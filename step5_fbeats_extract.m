function fbeats = step5_fbeats_extract(I2_orig,fs,f_r,gamma,alpha,win_len,win_type, lam_beat,HR_ma)
lam_curve = 20;
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
I2_orig_len = floor(length(I2_orig)/200)*200;
I2 = resample(I2_orig,100,fs);
basicTF.win = round(basicTF.win*3/5) ;

gg = I2; gg = gg-mean(gg);
%% apply the de-shape on the rough fECG to get the fetal HR
[~, ~, ~, tfrrF, ~, ~, tfrtic, t] = CFPH(gg, basicTF, advTF, cepR, P);

%% supress the possible residure of the mECG
for ti = 1:size(tfrrF, 2)
    idx = [round(HR_ma(ti)*0.94):round(HR_ma(ti)*1.06)] ;
    if max(HR_ma>=size(tfrrF,1))
        continue;
    else
        tfrrF(idx, ti) = tfrrF(idx, ti)/10 ;
    end
    
end

HR_fe = Extract_fhr(tfrrF, basicTF, lam_curve);
HR_fe3 = interp1(200:200:length(I2_orig), HR_fe, 1:length(I2_orig),'pchip','extrap') ;


%% apply the beat tracking to get the fetal HR
flocsf1p = beat_simple(I2, 100, HR_fe3.*basicTF.fr, lam_beat);
flocsf1q = beat_simple(-I2, 100, HR_fe3.*basicTF.fr, lam_beat);


%% get the fetal polarity
if abs(median(I2(flocsf1p))) > abs(median(I2(flocsf1q)))
    Po = 1 ; flocsf = flocsf1p ;
else
    Po = -1 ; flocsf = flocsf1q ;
    fprintf('\t\t*** reverse the fetal pole\n') ;
end

I2 = Po.*I2;
I2_orig = Po.*I2_orig ;

flocsf = RRconstraint(flocsf, I2, 100, 0.25);

I2 = I2_orig;
flocsf = round(flocsf*10/1);
fbeats = [];
SearchLen_p = round(48./mean(HR_fe3*basicTF.fr)) ;

for ii = 1:length(flocsf)
    [~, idx] = max(I2(max([1 flocsf(ii)-SearchLen_p]):min([flocsf(ii)+SearchLen_p length(I2)])));
    fbeats(ii) = flocsf(ii)-SearchLen_p+idx-1;
end
%fbeats = fbeats(fbeats>0); fbeats = fbeats(fbeats<length(I2_orig));
ind = find(fbeats>0 & fbeats<I2_orig_len);
fbeats = fbeats(ind);

I2_orig = Po.*I2_orig ;



end



