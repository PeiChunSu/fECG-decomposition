function Rpeaks = step2_Rpeaks(x,IHR,fs,fr,lam_beat)
basicTF.hop = 20;
basicTF.fs = 100;
basicTF.fr = fr;
basicTF.feat = 'SST11';
x0_len = length(x);
HR_ma2 = interp1(200:200:x0_len, IHR, 1:x0_len,'pchip','extrap') ;
%% Get maternal R peaks by DP
x1 = resample(x,basicTF.fs,fs);

mlocsf1p = beat_simple(x1, basicTF.fs, HR_ma2.*basicTF.fr, lam_beat);
mlocsf1q = beat_simple(-x1, basicTF.fs, HR_ma2.*basicTF.fr, lam_beat);

if abs(median(x1(mlocsf1p))) > abs(median(x1(mlocsf1q)))
    Po = 1 ; mlocsf1 = mlocsf1p ;
else
    Po = -1 ; mlocsf1 = mlocsf1q ;
end
x1 = Po.*x1 ;
x = Po.*x ;
mlocsf1 = RRconstraint(mlocsf1, x, 100, 0.25);

x_up = x;
mlocsf = round(mlocsf1*10);
Rpeaks = [];
SearchLen_p = round(48./mean(HR_ma2*basicTF.fr)) ;
SearchLen_q = 48;
for ii = 1:length(mlocsf)
    [~, idx] = max(x_up(max([mlocsf(ii)-SearchLen_p 1]):min([mlocsf(ii)+SearchLen_p length(x_up)]))) ;
    if ~isempty(idx)
        Rpeaks(ii) = mlocsf(ii)-SearchLen_p+idx-1 ;
    end
end
ind = find(Rpeaks>0 & Rpeaks<x0_len);
Rpeaks = Rpeaks(ind);

