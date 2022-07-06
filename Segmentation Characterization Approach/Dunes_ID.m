function [Skeleton_LC,GmphLC] = Dunes_ID(Gmph,MP,Conn,nhood,SK)

%Projet de Doctorat Willian Ney CASSOL
%Début de cette version du code 03 juin 2020

%Cette fonction vise la skeletonisation des lignes de crête des dunes à
%partir de:
    %Gmph = Geomorphons appliqué pour les lignes de crête
    %MP = %MP = Nombre minimale des pixels dans les objets - Lignes de crête
    %Conn = connectivité des pixels
    %nhood = Structuring element neighborhood
    %SK = filtre de skeletonisation

%Binarizer les geomorphons en considérant les classes Summit et Crest pour
%la ligne des crêtes des dunes
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

%Premier nettoyage des geomorphons binarizés
GmphLC = bwareaopen(GmphLC,MP,Conn);

%Effectuer une fermeture (Dilatation suivie d'une erosion)
GmphLC = imclose(GmphLC,nhood);
%Skeletonisation des zones de crêtes et sommet détectées par les
%geomorphons
Skeleton_LC = bwskel(GmphLC,'MinBranchLength',SK);

