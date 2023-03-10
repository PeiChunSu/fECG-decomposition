function  [Om0,Om0_real] = ECG_shrinkage0_color( x0,x0_real, current_beats, OptimalShrinkageOpt)

% optimal shrinakge of ECG

method = 'imp';
RRI = diff(current_beats);
RRI = [RRI(1) RRI];
MaximalQTp = ceil(quantile(RRI,0.95)*3/8);
MaximalQTt = ceil(quantile(RRI,0.95)*5/8);

tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;

tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;



V=[];
V2 = [];
II = [];
R = [];
S = [];
tRS = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    V2(:,count) = x0_real(idx);
    II = [II;idx];

end

[n_t, n_theta] = size(V);

cm = 1000; c= 4;
k = 20;
if n_theta>n_t
    if OptimalShrinkageOpt == "SN"
        [XN0, ~, r_p] = adaptiveHardThresholding(V, k, 'i');
    else
        [XN0,~,r_p] = optimal_shrinkage_color5(V,OptimalShrinkageOpt,'imp');
    end
else
    if OptimalShrinkageOpt == "SN"
        [XN0, ~, r_p] = adaptiveHardThresholding(V', k, 'i');
    else
        [XN0,~,r_p] = optimal_shrinkage_color5(V',OptimalShrinkageOpt,'imp');
    end
    XN0 = XN0';
end

if n_theta>n_t
    if OptimalShrinkageOpt == "SN"
        [XN02, ~, r_p] = adaptiveHardThresholding(V2, k, 'i');
    else
        [XN02,~,r_p] = optimal_shrinkage_color5(V2,OptimalShrinkageOpt,'imp');
    end
else
    if OptimalShrinkageOpt == "SN"
        [XN02, ~, r_p] = adaptiveHardThresholding(V2', k, 'i');
    else
        [XN02,~,r_p] = optimal_shrinkage_color5(V2',OptimalShrinkageOpt,'imp');
    end
    XN02 = XN02';
end








Om0 = zeros(1,length(x0));
Om0_real = zeros(1,length(x0_real));

for s = 1:length(current_beats)

    XN_to_add = XN0(:,s);
    XN2_to_add = XN02(:,s);
    %reconstruction

    if s == 1
        left_overlap = 2;
        right_overlap = length(intersect(II(s+1,:),II(s,:)));
    elseif s == length(current_beats)
        left_overlap = length(intersect(II(s-1,:),II(s,:)));
        right_overlap = 2;
    else
        left_overlap = length(intersect(II(s-1,:),II(s,:)));
        right_overlap = length(intersect(II(s+1,:),II(s,:)));
    end


    if left_overlap <= 1; left_overlap = 2; end
    if right_overlap <= 1; right_overlap = 2; end


    W = ones(MaximalQTp+MaximalQTt+1,1);
    W(1:left_overlap) = sin(linspace(0,pi/2,left_overlap)).^2;
    W(end:-1:end-right_overlap+1) = sin(linspace(0,pi/2,right_overlap)).^2;
    %Om_toadd = interp1(II(s,:),XN(:,loc),II(s,:));
    %Om_toadd = W'.*Om_toadd;
    Om_toadd = (W.*XN_to_add)';
    Om_real_toadd = (W.*XN2_to_add)';

    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
    Om0_real(II(s,:)) = Om0_real(II(s,:)) + Om_real_toadd;

end

Om0(Om0 == 0) = x0(Om0 == 0);
Om0_real(Om0_real == 0) = x0_real(Om0_real == 0);







