function [pente] =pente_app(Dune_Profil)
d2r = pi/180;
%Calcul de la différence entre le fond plat et le plan d'inclinaison pente°
[m,n] = size(Dune_Profil);
Y = Dune_Profil(:,1);
Z = Dune_Profil(:,2);

for i=2:m
    pente(i) = atand((Z(i)-Z(i-1))/(Y(i)-Y(i-1)));
 end   



