function [Om,Of,mbeats,fbeats] = step6_iteration(x,Of0,fs,fr,gamma,alpha,win_len,win_type,lam_beat,OptimalShrinkageOpt)
x_rm = x-Of0;
tfrr = step1_deshape(x_rm,fs,fr,gamma,alpha,win_len,win_type);
IHR = step1_IHR(tfrr,fr);
mbeats = step2_Rpeaks(x_rm,IHR,fs,fr,lam_beat);
XN0 = step3_OS( x_rm,mbeats, OptimalShrinkageOpt);
Om = step4_stitch(XN0,x_rm,mbeats);
[Of,fbeats] = step5_fECG(x,Om,fs,fr,gamma,alpha,win_len,win_type,lam_beat,OptimalShrinkageOpt,IHR);
