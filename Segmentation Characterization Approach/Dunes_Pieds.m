function [GmphBD,BDSelec,Skelette_LC_Final,C,L] = Dunes_Pieds(GmphBD0,DTM,Skeleton_LC,MP,Conn)

%Projet de Doctorat Willian Ney CASSOL
%D�but de cette version du code 03 juin 2020

%Cette fonction vise la skeletonisation des lignes de cr�te des dunes �
%partir de:


Vector_gmph = [1 2 3 4 5 6 7 8 9 10];
%Binarizer les geomorphons en consid�rant les classes 8=footslope 9=valley
%10=depression 1=flat 7=hollow
[k,l] = size(GmphBD0);
for i=1:k
    for j=1:l
    if (GmphBD0(i,j)==Vector_gmph(9)||GmphBD0(i,j)==Vector_gmph(10)||GmphBD0(i,j)==Vector_gmph(1)||GmphBD0(i,j)==Vector_gmph(7)||GmphBD0(i,j)==Vector_gmph(8))
        GmphBD(i,j)=1;        
    else
        GmphBD(i,j)=0;   
    end
    end
end

%Premier nettoyage des geomorphons binariz�s pour les pieds des dunes
GmphBD = bwareaopen(GmphBD,MP);

%Trouver le p�rim�tre des dunes
BDSelec = bwperim(GmphBD,8);
GmphBD = imcomplement(GmphBD);
[k,l] = size(GmphBD);
for i=1:k
    for j=1:l
    if (DTM(i,j)== 0)
       GmphBD(i,j)= 0;   
    else
    end
    end
end


%Filtrer le skelete des lignes de cr�te pour seulement celles se trouvant
%dans les r�gions des dunes - Pieds des dunes
[k,l] = size(GmphBD);
for i=1:k
    for j=1:l
    if (GmphBD(i,j)==1)
        Skelette_LC_Final(i,j)=Skeleton_LC(i,j);
    else
        Skelette_LC_Final(i,j)=0;   
    end
    end
end

[k,l] = size(Skelette_LC_Final);
for i=1:k
    for j=1:l
    if (DTM(i,j)== 0)
       Skelette_LC_Final(i,j)= 0;   
    else
    end
    end
end


%Chercher les coordonn�es des lignes skeletonis�es afin de les utiliser
%pour ensuite conserver les objets dunes dans leur d�tection - Cela permet
%de garder seulement les "objets dunes" ayant une ligne de cr�te dans la
%r�gion des dunes d�limit�es par les pieds des dunes
[k,l] = size(Skelette_LC_Final);
w = 1;
for i=1:k
    for j=1:l
    if (Skelette_LC_Final(i,j)==1)
        Skelette_LC_Final(i,j)=Skeleton_LC(i,j);
        C(w) = j;
        L(w) = i;
        w=w+1;
    else
    end
    end
end

%Conserver seulement les "objets dunes" ayant une ligne de cr�te en
%superposition
BDSelec = bwselect(GmphBD,C,L,8);
[k,l] = size(BDSelec);
for i=1:k
    for j=1:l
    if (DTM(i,j)== 0)
       BDSelec(i,j)= 0;   
    else
    end
    end
end
