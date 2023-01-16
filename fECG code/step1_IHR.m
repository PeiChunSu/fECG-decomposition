function IHR = step1_IHR(tfrr,f_r)
lam_curve = 20;
basicTF.hop = 20;
basicTF.fs = 100;
basicTF.fr = f_r;
basicTF.feat = 'SST11';
IHR = Extract_mhr(tfrr, basicTF, lam_curve);