function Om0 = step4_stitch(XN0,x0,current_beats)

RRI = diff(current_beats);
RRI = [RRI(1) RRI];
MaximalQTp = ceil(quantile(RRI,0.95)*3/8);
MaximalQTt = ceil(quantile(RRI,0.95)*5/8);

tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;

II = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    II = [II;idx];
end
Om0 = zeros(1,length(x0));

for s = 1:length(current_beats)

    XN_to_add = XN0(:,s);
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
    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
end

Om0(Om0 == 0) = x0(Om0 == 0);



