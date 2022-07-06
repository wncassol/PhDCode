function [Length,Dc] = Ditances(B_DLC_Break,ID_Coord_All)


%calculer distance geodesique
%ID_Coord_All = ID Ligne Crete, Xd, Yd, Xf, Yf
ID = ID_Coord_All(:,1);
Xd = ID_Coord_All(:,2);
Yd = ID_Coord_All(:,3);
Xf = ID_Coord_All(:,4);
Yf = ID_Coord_All(:,5);
%Geodetic distance
Dc = sqrt((Xf-Xd).^2+(Yf-Yd).^2);
    
%Legth estimatiojn 
n_obj = length(B_DLC_Break);
for i=1:n_obj
    C = B_DLC_Break{i,:};
    [Xo] = C(:,1);
    [Yo] = C(:,2);
    pos =1;
    w = length(Xo);
    for j=2:w
        dist(j) = sqrt((Xo(j)-Xo(j-1))^2+((Yo(j)-Yo(j-1))^2));
    end
    Dist(i) = sum(dist)/2;
    dist=0;
    
end
% CrestlineLength
Length = Dist;
    
