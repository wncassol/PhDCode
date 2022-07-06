function [Skeleton_LC,GmphLC] = Dunes_ID(Gmph,MP,Conn,nhood,SK)

%Projet de Doctorat Willian Ney CASSOL
%D�but de cette version du code 03 juin 2020

%Cette fonction vise la skeletonisation des lignes de cr�te des dunes �
%partir de:
    %Gmph = Geomorphons appliqu� pour les lignes de cr�te
    %MP = %MP = Nombre minimale des pixels dans les objets - Lignes de cr�te
    %Conn = connectivit� des pixels
    %nhood = Structuring element neighborhood
    %SK = filtre de skeletonisation

%Binarizer les geomorphons en consid�rant les classes Summit et Crest pour
%la ligne des cr�tes des dunes
Vector_gmph = [1 2 3 4 5 6 7 8 9 10];
[k,l] = size(Gmph);
for i=1:k
    for j=1:l
    if (Gmph(i,j)== Vector_gmph(2) || Gmph(i,j)== Vector_gmph(3))
        GmphLC(i,j)=1;        
    else
        GmphLC(i,j)=0;   
    end
    end
end

%Premier nettoyage des geomorphons binariz�s
GmphLC = bwareaopen(GmphLC,MP,Conn);

%Effectuer une fermeture (Dilatation suivie d'une erosion)
GmphLC = imclose(GmphLC,nhood);
%Skeletonisation des zones de cr�tes et sommet d�tect�es par les
%geomorphons
Skeleton_LC = bwskel(GmphLC,'MinBranchLength',SK);

