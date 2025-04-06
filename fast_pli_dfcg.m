function dfcg=fast_pli_dfcg(multi,window1,step1)


%%% A FAST IMPLEMENTATION OF PLI FOR DYNAMIC FUNCTIONAL BRAIN NETWORKS
%%% INPUT : multi = filtered multichannel recordings with dimensions 
%%%                 sensors/sources/rois x samples (time points)
%%%       window1 = width of temporal segment in samples
%%%         step1 = moving step of the center of temporal window in samples
%%% OUTPUT : dfcg = temporal segments x sensors/sources/rois x sensors/sources/rois

%STAVROS I. DIMITRIADIS 17/05/2018
% CARDIFF UNIVERSITY BRAIN RESEARCH IMAGING CENTRE (CUBRIC)
% Neuroinformatics Group, CUBRIC, CARDIFF,WALES,UK
%http://users.auth.gr/~stdimitr/index.html

[rois samples]=size(multi);


slides=round((samples-window1)/step1);

phases=angle(hilbert(multi'))';



dfcg=zeros(slides,rois,rois);

for k=1:rois
    for l=(k+1):rois
        df=wrapToPi(phases(k,:) - phases(l,:));
            for ts=1:slides
                tt=[(ts-1)*step1 + 1:(ts-1)*step1 + window1];
                pli=abs(mean(sign(df(tt))));
                dfcg(ts,k,l)=pli;
            end
    end
end



