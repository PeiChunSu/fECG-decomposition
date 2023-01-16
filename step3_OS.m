function  XN0 = step3_OS( x0,current_beats, OptimalShrinkageOpt)

% optimal shrinakge of ECG


RRI = diff(current_beats);
RRI = [RRI(1) RRI];
MaximalQTp = ceil(quantile(RRI,0.95)*3/8);
MaximalQTt = ceil(quantile(RRI,0.95)*5/8);

tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;


V=[];
II = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    II = [II;idx];

end

[n_t, n_theta] = size(V);

if n_theta>n_t
    [XN0,~,r_p] = optimal_shrinkage_color5(V,OptimalShrinkageOpt,'imp');

else

    [XN0,~,r_p] = optimal_shrinkage_color5(V',OptimalShrinkageOpt,'imp');

    XN0 = XN0';
end

