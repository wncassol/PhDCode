function [OBD_Plat,OBD_NonPlat,OBD_NonPlat_Skeleton,Limites_Bas_Dunes] = NonDunes_Objects(L_OND,gmph_out,NonDunes_Clean2,SK)
%Projet de Doctorat Willian Ney CASSOL
%Début de cette version du code 02 juin 2020

%Fonction pour le retour des limites des pieds des dunes
%Les entrants:
    %LWTF - Image avec les différents objets identifiés dans les dunes avec la
    %fonction [BWTF,LWTF,nWTF,AWTF] =  bwboundaries(NDSelec,Conn);
    %gmph_out - Geomorphons généré à partir des régions des bas de dunes
    %dientifiées
    %NonDunes_Clean2 - Region des pieds des dunes identifiés et
    %binarizés
    %SK -  Minimum branch length
%Les sortants:
    %OBD_Plat - Image avec les objets non-dunes considérés plats
    %OBD_NonPlat - Image avec les objets non-dunes considérés non-plats
    %OBD_NonPlat_Skeleton - Skeleton des OBD_NonPlat en considérant la
    %fonction bwskel et SK
    %Limites_Bas_Dunes - Image avec la fusion des OBD_Plat et OBD_NonPlat_Skeleton

%Stockage dans un vecteur du nombre d'objets ''non-dune'' de la zone
Vector = [0 ];
Wi = 2;
[k,l] = size(L_OND);
for i=1:k
    for j=1:l
    if (L_OND(i,j)~=0)
        if (L_OND(i,j)~=Vector(Wi-1))
            Vector(Wi)=L_OND(i,j); 
            Wi = Wi + 1;
        else
        end
    else
    end
    end
end
Vector= unique(Vector,'sort');
Vector(1)=[];
%Nombre total des pixels de chaque "Objet non-dune trouvé"
[kW,lW] = size(Vector);
for i=1:lW
    Nombre_Pixels_ObNonDune(i) = sum(L_OND(:)==Vector(i));
end

Vector_Nombre_Pixels_ObNonDune(1,:) = [Vector];
Vector_Nombre_Pixels_ObNonDune(2,:) = [Nombre_Pixels_ObNonDune];
Vector_Nombre_Pixels_ObNonDune(3,:) = zeros(size(Vector_Nombre_Pixels_ObNonDune(2,:)));
Vector_Nombre_Pixels_ObNonDune(4,:) = zeros(size(Vector_Nombre_Pixels_ObNonDune(2,:)));
% Vector_Nombre_Pixels_ObNonDune(5,:) = zeros(size(Vector_Nombre_Pixels_ObNonDune(2,:)));

%Les dix geomorphons possibles
Vector_gmph = [1 2 3 4 5 6 7 8 9 10];
[k,l] = size(L_OND);
[kW,lW] = size(Vector);
for iW=1:lW
for i=1:k
    for j=1:l
        if (L_OND(i,j)==Vector(iW))
            if (gmph_out(i,j)==Vector_gmph(7) || gmph_out(i,j)==Vector_gmph(9) || gmph_out(i,j)==Vector_gmph(10))% gmph_out(i,j)==Vector_gmph(8) || gmph_out(i,j)==Vector_gmph(6))
                
                Vector_Nombre_Pixels_ObNonDune(3,iW) = Vector_Nombre_Pixels_ObNonDune(3,iW)+1;
                
            elseif (gmph_out(i,j)==Vector_gmph(1))% ||  gmph_out(i,j)==Vector_gmph(6)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Vector_Nombre_Pixels_ObNonDune(4,iW) = Vector_Nombre_Pixels_ObNonDune(4,iW)+1;

            else
%                 Vector_Nombre_Pixels_ObNonDune(5,iW) = Vector_Nombre_Pixels_ObNonDune(5,iW)+1;
            end
        else
        end
    end
end
end

%Détermination si l'objet non-dune est plat ou non
Vector_Nombre_Pixels_ObNonDune(5,:) = zeros(size(Vector_Nombre_Pixels_ObNonDune(2,:)));
[kW,lW] = size(Vector_Nombre_Pixels_ObNonDune);
for i=1:lW
    if(Vector_Nombre_Pixels_ObNonDune(3,i)>Vector_Nombre_Pixels_ObNonDune(4,i))
        % 1= Plat / 0 = Non plat
        Vector_Nombre_Pixels_ObNonDune(5,i) = 0;
    else
         % 1= Plat / 0 = Non plat
        Vector_Nombre_Pixels_ObNonDune(5,i) = 1;
    end
end
%En fonction de la détermination précedente, on détermine que si l'objet
%non-dune est plat, il n'y a pas de skeletonisation et la frontière de
%cette région est utilisée comme pied de dune. Sinon, si elle n'est pas
%plat (non plat), elle sera skeletonisée et cette skelete sera considéré
%comme la frontière de la dune, ou pied de dune.

%Création de deux nouvelles images, une à sekeletoniser et l'autre à
%conserver.
[k,l] = size(L_OND);
[kW,lW] = size(Vector);
OBD_NonPlat = zeros(size(gmph_out));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iW=1:lW
for i=1:k
    for j=1:l
        if (L_OND(i,j)==Vector(iW))
                     % 1= Plat / 0 = Non plat
            if (Vector_Nombre_Pixels_ObNonDune(5,iW)==0) 
               OBD_NonPlat(i,j)=NonDunes_Clean2(i,j);
            else
               OBD_NonPlat(i,j)=0;

            end
        else
        end
    end
end
end
%Les objets non-dunes qui ne sont pas plats, doivent être skeletonisés pour
%trouver les lignes qui délimitent les pieds des dunes:
OBD_NonPlat_Skeleton_bin = imbinarize(OBD_NonPlat);
OBD_NonPlat_Skeleton = bwskel(OBD_NonPlat_Skeleton_bin,'MinBranchLength',SK);

%Les objets non-dunes considérés comme plats, ne doivent pas être 
%skeletonisés pour trouver les lignes qui délimitent les pieds des dunes,
%au contraire, l'objet même et ses frontières représentent une zone
%non-dune qui n'appartient à aucune structure et ainsi les limites
%(frontières) de cet objet sont considérées comme les bas des dunes
%délimitées par les lignes de crêtes
[k,l] = size(L_OND);
[kW,lW] = size(Vector);
OBD_Plat = zeros(size(gmph_out));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iW=1:lW
for i=1:k
    for j=1:l
        if (L_OND(i,j)==Vector(iW))
                     % 1= Plat / 0 = Non plat
            if (Vector_Nombre_Pixels_ObNonDune(5,iW)==1) 
               OBD_Plat(i,j)=NonDunes_Clean2(i,j);
            else
               OBD_Plat(i,j)=0;

            end
        else
        end
    end
end
end

% Limites_Bas_Dunes = imfuse(OBD_NonPlat_Skeleton,OBD_Plat,'falsecolor','Scaling','joint');
[k,l] = size(L_OND);
for i=1:k
    for j=1:l
        if (OBD_Plat(i,j)==1 || OBD_NonPlat_Skeleton(i,j)==1)
            Limites_Bas_Dunes(i,j)=1;%%%%%%%%%%%%%%%%
        else
            Limites_Bas_Dunes(i,j)=0;
        end
    end
end