function D = pearson_cc(X,Y) % calculation of distance ...
        
%%% INPUT: X,Y variables
%%% OUTUT: pearson correlation value

%STAVROS I. DIMITRIADIS 28/01/2012

%% Contact: DimitriadisS@cardiff.ac.uk/ stidimitriadis@gmail.com
%% WEBPAGE: https://www.researchgate.net/profile/Stavros_Dimitriadis

        X=X';
        Y=Y';
        N=length(X);
        XY=X.*Y;
        X2 = X.*X;
        Y2=Y.*Y;
        num = (N.*(sum(XY)))-(sum(X).*sum(Y));
        densq = (N.*sum(X2)-(sum(X).^2)).*(N.*sum(Y2)-(sum(Y).^2));
        den = sqrt(densq);
        D = num./den;