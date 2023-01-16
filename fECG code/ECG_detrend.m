function x0 = ECG_detrend(x0, num_med_s, num_med_l)
%detrend
%if_morph = 1 : preserve morphology
%         = 0 : otherwise
if nargin < 5
    if_smooth = 1;
end
    

x0_trend = movmedian(x0,num_med_s);% medfilt2(x0, [1, num_med_s]);
x0_trend = movmedian(x0_trend,num_med_l);% medfilt2(x0_trend, [1, num_med_l]);
x0 = x0 - x0_trend;
