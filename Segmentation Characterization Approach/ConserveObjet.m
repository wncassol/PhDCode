function [BDSelec] = ConserveObjet(Skelette_LC_Final,Dunes_Clean1)

%Projet de Doctorat Willian Ney CASSOL
%D�but de cette version du code 03 juin 2020
%Chercher les coordonn�es des lignes skeletonis�es afin de les utiliser
%pour ensuite conserver les objets dunes dans leur d�tection - Cela permet
%de garder seulement les "objets dunes" ayant une ligne de cr�te dans la
%r�gion des dunes d�limit�es par les pieds des dunes
[k,l] = size(Skelette_LC_Final);
w = 1;
for i=1:k
    for j=1:l
    if (Skelette_LC_Final(i,j)==1)
        C(w) = j;
        L(w) = i;
        w=w+1;
    else
    end
    end
end

%Conserver seulement les "objets dunes" ayant une ligne de cr�te en
%superposition
BDSelec = bwselect(Dunes_Clean1,C,L,8);
[k,l] = size(Dunes_Clean1);
