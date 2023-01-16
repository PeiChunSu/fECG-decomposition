function [Of0,fbeats] = step5_fECG(x,Om0,fs,fr,gamma,alpha,win_len,win_type,lam_beat,OptimalShrinkageOpt,HR_ma)

I2_orig = x-Om0;
fbeats = step5_fbeats_extract(I2_orig,fs,fr,gamma,alpha,win_len,win_type, lam_beat,HR_ma);
XN0 = step3_OS(I2_orig,fbeats, OptimalShrinkageOpt);
Of0 = step4_stitch(XN0,I2_orig,fbeats);