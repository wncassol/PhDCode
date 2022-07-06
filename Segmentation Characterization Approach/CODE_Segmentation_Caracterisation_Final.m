%Approche de segmentation et de caractérisation des dunes sous-marines
%Implémentation de l'approche développée dans le projet de Doctorat de
%Willian Ney CASSOL
%Version finale 06 juillet 2022
%La surface de test (G15_19_1_100.tif) a été fournie par groupe océan et le service Hydrographique du Canada. Ces données doivent être seulement utilisées pour la validation de ce code.
%La surface G15_19_1_100_gmph14.tif a été générée avec l'algorithme des
%geomorphons par la fonction 
clear all
close all
format compact
format long g
clc
tic
%--------------------------------------------------------------------------
%Paramètres utilisés dans les traitements pour la délimitation des lignes
%de crête et les pieds des dunes à partir des Geomorphons.
%SK = filtre de skeletonisation - Taille minimale de la branche
%squeletonisée conservée
SK = 10;
%Conn = connectivité des pixels
Conn = 8;
%MP = Nombre minimale des pixels dans les objets - Lignes de crête
MP1 = 10;
%MP = Nombre minimale des pixels dans les objets dunes
MP = 50;
%nhood = Structuring element neighborhood
nhood = [1 1 1; 1 1 1; 1 1 1];
%Lecture du Modèle Numérique de Terrain - Lecture pour la matrice de
%géoréférencement DTM *****************************************************
[DTM,R1DTM,R2DTM] = geotiffread('G15_19_1_100.tif');
DTM(DTM==-9999)=0;

%Lecture des surfaces geomorphon ******************************************
GmphLC0 = geotiffread('G15_19_1_100_gmph14.tif');
GmphBD0 = GmphLC0;

%Identification des lignes de crête et skeletonisation à partir des
%Geomorphons GmphLC0 et de la fonction Dunes_ID
%Selon le diagramme avec le modële opêrationnel - 2ème étape
[Skeleton_LC,GmphLC] = Dunes_ID(GmphLC0,10,Conn,nhood,SK);
%%%Astuce pour enlever les "loops" dans les lignes de crête
Skeleton_LC = imfill(Skeleton_LC,'holes');
Skeleton_LC = bwskel(Skeleton_LC);

%Identification des bas des dunes et ainsi les objets dunes primodiales à
%l'aide de la fonction Dunes_Pieds
%Selon le diagramme avec le modële opêrationnel - 3ème étape
[GmphBD,BDSelec,Skelette_LC_Final,C,L] = Dunes_Pieds(GmphBD0,DTM,Skeleton_LC,MP,Conn);

%Selon le diagramme avec le modële opêrationnel - 4ème étape
%Ouverture morphologique de l'image avec les éléments des dunes
%primordiales
Dunes_primordiales = imopen(BDSelec,nhood);
% figure;imshow(Dunes_primordiales);title('Ouverture morphologique')

%--------------------------------------------------------------------------
%Début des opérations avec les objets dunes

%Première application des contours actifs pour les objets dunes
%primordiales

%Filtrer le skeleton des lignes de crête pour seulement celles se trouvant
%dans les régions des dunes - Pieds des dunes - après la première itération
%des contours actifs
[k,l] = size(Dunes_primordiales);
for i=1:k
    for j=1:l
    if (Dunes_primordiales(i,j)==1)
        Skelette_LC_Final(i,j)=Skelette_LC_Final(i,j);
    else
        Skelette_LC_Final(i,j)=0;   
    end
    end
end

NonDunes = imcomplement(Dunes_primordiales);
Dunes_primordiales = bwareaopen(Dunes_primordiales,MP,Conn);
NonDunes = imcomplement(Dunes_primordiales);

[k,l] = size(Dunes_primordiales);
for i=1:k
    for j=1:l
    if(DTM(i,j)==0)
        NonDunes(i,j)=0; 
    else
    end
    end
end

%Nettoyage des petits "objets"(pixels) qui ne sont pas des entre dunes
Non_Dunes_Clean = bwareaopen(NonDunes,MP1,Conn);
[k,l] = size(Non_Dunes_Clean);
for i=1:k
    for j=1:l
    if(DTM(i,j)==0)
        Non_Dunes_Clean(i,j)=0; 
    else
    end
    end
end

%Skeletonisation des zones de vallée détectées par les geomorphons
Skeleton_Limites_Dunes = bwskel(Non_Dunes_Clean,'MinBranchLength',SK);

[k,l] = size(Skeleton_Limites_Dunes);
w = 1;
for i=1:k
    for j=1:l
    if (Skeleton_Limites_Dunes(i,j)==1)
        C(w) = j;
        L(w) = i;
        w=w+1;
    else
    end
    end
end

%Conserver seulement les "objets non-dunes" ayant une ligne de bas de dunes en
%superposition
Objets_Non_Dunes_Selec = Non_Dunes_Clean;

%--------------------------------------------------------------------------
%Début du traitement des images avec des objets dunes avec la fonction bwboundaries
%Légende: OND
[B_OND,L_OND,n_OND,A_OND] =  bwboundaries(Objets_Non_Dunes_Selec,Conn);

[BTEST,LTEST] = bwboundaries(Objets_Non_Dunes_Selec,Conn,'noholes');

L_OND = LTEST;

[k,l] = size(L_OND);
for i=1:k
    for j=1:l
    if (L_OND(i,j)>=1)
       NonDunes_Clean1(i,j)=1; 
    else
        NonDunes_Clean1(i,j)=0; 
    end
    end
end

%Supprimer tous les objets plus petits que MP1 pixels de l'image
NonDunes_Clean1 = bwareaopen(L_OND,MP1,Conn);
[k,l] = size(L_OND);
for i=1:k
    for j=1:l
    if (DTM(i,j)==0)
       NonDunes_Clean1(i,j)=0;
           else
    end
    end
end

Dunes_Clean1 = imcomplement(NonDunes_Clean1);
[k,l] = size(L_OND);
for i=1:k
    for j=1:l
    if (DTM(i,j)==0)
       Dunes_Clean1(i,j)=0;
           else
    end
    end
end

Dunes_Clean1 = imopen(Dunes_Clean1,nhood);
%Filtrer squelettes des dunes mineurs que MP1
Skelette_LC_Final = bwareaopen(Skelette_LC_Final,MP1,Conn);
%Filtrer le skeleton des lignes de crête pour seulement celles se trouvant
%dans les régions des dunes - Pieds des dunes
[k,l] = size(Dunes_Clean1);
for i=1:k
    for j=1:l
    if (Dunes_Clean1(i,j)==1)
        Skelette_LC_Final(i,j)=Skelette_LC_Final(i,j);
    else
        Skelette_LC_Final(i,j)=0;   
    end
    end
end
Dunes_Clean1 = ConserveObjet(Skelette_LC_Final,Dunes_Clean1);

%Selon le diagramme avec le modële opêrationnel - 5ème étape et 6ème étape
NonDunes_Clean2 = imcomplement(Dunes_Clean1);

[k,l] = size(NonDunes_Clean2);
for i=1:k
    for j=1:l
    if (DTM(i,j)==0)
       NonDunes_Clean2(i,j)=0;
           else
    end
    end
end

[B_OND,L_OND,n_OND,A_OND] =  bwboundaries(NonDunes_Clean2,Conn);
%OND = Objet non Dune
[BTEST,LTEST] = bwboundaries(NonDunes_Clean2,Conn,'noholes');
L_OND = LTEST; 

%Nettoyage des petits objets avec une aire inférieure à MP
Dunes_Clean1 = bwareaopen(Dunes_Clean1,MP,Conn);
NonDunes_Clean1 = bwareaopen(NonDunes_Clean1,MP,Conn);
NonDunes_Clean2 = bwareaopen(NonDunes_Clean2,MP,Conn);

%%%%%%%DUNES OBJECTS ------------------------------------------------------
[OBD_Plat,OBD_NonPlat,OBD_NonPlat_Skeleton,Limites_Bas_Dunes] = NonDunes_Objects(L_OND,GmphBD0,NonDunes_Clean2,SK);

Limites_Bas_Dunes = imcomplement(Limites_Bas_Dunes);
[k,l] = size(Limites_Bas_Dunes);
for i=1:k
    for j=1:l
    if (DTM(i,j)==0)
       Limites_Bas_Dunes(i,j)=0;
           else
    end
    end
end

Dunes_Clean2 = imfill(Dunes_Clean1,'holes');

% Enregistrement des images avec les limites estimês des dunes
coordRefSysCode = 32187;
%Image avec les pieds des dunes
FINAL_PIEDSDUNES = Limites_Bas_Dunes;
FINAL_PIEDSDUNES = imbinarize(FINAL_PIEDSDUNES);

%Image avec la ligne de crete
FINAL_LIGNESCRETE = Skelette_LC_Final;
%%%Astuce pour enlever les "loops" dans les lignes de crête
FINAL_LIGNESCRETE = imfill(FINAL_LIGNESCRETE,'holes');
FINAL_LIGNESCRETE = bwskel(FINAL_LIGNESCRETE);

%***************************************************************************************************************************************
geotiffwrite('G15_19_1_100_LigneCrete.tif',FINAL_LIGNESCRETE, R1DTM,'CoordRefSysCode', coordRefSysCode)


%Génération seulement des lignes limites des dunes et non des objets et
%polygones
PERIM_BD_1 = bwperim(OBD_Plat,Conn);
[k,l] = size(PERIM_BD_1);
for i=1:k
    for j=1:l
    if (PERIM_BD_1(i,j)==1 || OBD_NonPlat_Skeleton(i,j)==1)
       PIEDSDUNES_PERIM(i,j)=1;
    else
        PIEDSDUNES_PERIM(i,j)=0;
    end
    end
end

%Image avec le perimètre du pied de dune
PIEDSDUNES_PERIM = imbinarize(PIEDSDUNES_PERIM);

[B_DPP,L_DPP,n_DPP,A_DPP] =  bwboundaries(PIEDSDUNES_PERIM,Conn);
n_dpp = length(B_DPP);
%Nettoyage  des pieds des dunes situés dans le contour de l'image binarisée
for i=1:n_dpp
    C = B_DPP{i};
    [Xdpp] = C(:,1);
    [Ydpp] = C(:,2);
    [k,l] = size(L_DPP);
    for j=1:length(Xdpp)
        x0 = Xdpp(j); 
        x1 = Xdpp(j)+1;    
        x2 = Xdpp(j)-1; 
        y0 = Ydpp(j);
        y1 = Ydpp(j)+1; 
        y2 = Ydpp(j)-1; 
        if ((x1>k) || (x2>k) || (y1>l) || (y2>l) || (x1==0) || (x2==0) || (y1==0) || (y2==0))
        else    
            if (((DTM(x2,y2))==0) || ((DTM(x2,y0))==0) || ((DTM(x2,y1))==0) || ((DTM(x0,y2))==0) || ((DTM(x0,y1))==0) || ((DTM(x1,y2))==0) || ((DTM(x1,y0))==0) || ((DTM(x1,y1))==0))
                PIEDSDUNES_PERIM(x0,y0)=0;
            else
            end
        end
    end
end

%***************************************************************************************************************************************
geotiffwrite('G15_19_1_100_PiedsDunesPerim.tif',PIEDSDUNES_PERIM, R1DTM,'CoordRefSysCode', coordRefSysCode)
figure;hold on; title('Pieds des dunes et Ligne de crête');imshow(imfuse(FINAL_LIGNESCRETE,PIEDSDUNES_PERIM));


%PHASE II  *****************************************************************************************************************
%--------------------------------------------------------------------------
%Paramètres utilisés dans les traitements pour la délimitation des lignes
%de crête et les pieds des dunes à partir des Geomorphons.
%SK = filtre de skeletonisation - Taille minimale de la branche
%squeletonisée conservée
SK = 10;
%Conn = connectivité des pixels
Conn = 8;
%MP = Nombre minimale des pixels dans les objets - Lignes de crête
MP1 = 20;
%MP = Nombre minimale des pixels dans les objets dunes
MP = 60;
%nhood = Structuring element neighborhood
nhood = ones(3);

%DPP = Dune_Pied_Perim
DPP = PIEDSDUNES_PERIM;
%DLC = Dune_Ligne_Crête
DLC = FINAL_LIGNESCRETE;
[B_DLC,L_DLC,n_DLC,A_DLC] =  bwboundaries(DLC,Conn);

%Selon le diagramme avec le modële opêrationnel - 1ère étape
branchimage = bwmorph(DLC, 'branchpoints');
DLC_Break = DLC-branchimage;
DLC_Break = bwareaopen(DLC_Break,MP1);
DLC_Break1 = DLC_Break+branchimage;
branchimage1 = bwmorph(DLC_Break1, 'branchpoints');
DLC_Break2 = DLC_Break1-branchimage1;
DLC_Break = DLC_Break2;
DLC_Break = bwareaopen(DLC_Break,MP1);
[B_DLC_Break,L_DLC_Break,n_DLC_Break,A_DLC_Break] =  bwboundaries(DLC_Break,Conn);

Debut_Fin_DLC = bwmorph(DLC_Break,'endpoints');
[row,col] = find(Debut_Fin_DLC);
for i=1:length(row)
        ID_Obj_DLC = L_DLC_Break(row(i),col(i));
        %coordonnées début et fin avec numéro de ligne de crête
        ID_Coord(i,:) = [ID_Obj_DLC,col(i),row(i)];
end

%Association de la coordonnée de début et fin d'une ligne de crête dans la
%même ligne de la matrice correspondante à l'objet
%Essaie de création du vecteur ID_Obj Xi Yi Xf Yf
n_obj = max(ID_Coord(:,1));
for i=1:n_obj
   %trouver le numéro de la ligne dans la matrice ID_Coord qui correspond
   %au début et à la fin de la ligne de crête
   Debut_fin_DLC(:,i) =  find(ID_Coord(:,1)==i);
end

ID_Coord_All = zeros(n_obj,5);
%Association des coordonnées de début et fin avec le numéro de la ligne de
%crÊte dans une seule matrice
for i=1:n_obj
    %Numéro de la ligne de crête
    ID_Coord_All(i,1) =  i;
    %Coordonnée X de début
    ID_Coord_All(i,2) =  ID_Coord(Debut_fin_DLC(1,i),2);
    %Coordonnée Y de début
    ID_Coord_All(i,3) =  ID_Coord(Debut_fin_DLC(1,i),3);
    %Coordonnée X de fin
    ID_Coord_All(i,4) =  ID_Coord(Debut_fin_DLC(2,i),2);
    %Coordonnée Y de fin
    ID_Coord_All(i,5) =  ID_Coord(Debut_fin_DLC(2,i),3);
end

%Imprimer les coordonnées avec les numéros de début et fin de crête
figure;
imshow(imfuse(DLC_Break,DPP));
%Nettoyage
DPP = bwareaopen(DPP,MP1,Conn);

hold on;plot(ID_Coord_All(:,2),ID_Coord_All(:,3),'r*');hold on
plot(ID_Coord_All(:,4),ID_Coord_All(:,5),'b*')
hold on;

Text = string(ID_Coord_All(:,1));
for i=1:length(ID_Coord_All(:,1))
text(ID_Coord_All(i,2),ID_Coord_All(i,3),Text(i),'Color','r','FontSize',14)
end
hold on
for i=1:length(ID_Coord_All(:,1))
text(ID_Coord_All(i,4),ID_Coord_All(i,5),Text(i),'Color','b','FontSize',14)
end
title('Lignes de crête avec les points de début et de fin'); 
hold on
for i=1:length(ID_Coord_All(:,1))
plot([ID_Coord_All(i,2),ID_Coord_All(i,4)],[ID_Coord_All(i,3),ID_Coord_All(i,5)],'Color','green')
end

%Suppréssion des objets inférieurs à MP1
DPP = bwareaopen(DPP,MP1,Conn);

%rélocalisation des objets dans l'image sans les pieds des dunes dans le
%bord de l'image et sans objets inférieurs à MP1
[B_DPP,L_DPP,n_DPP,A_DPP] =  bwboundaries(DPP,Conn);

%éléments structurants au cas ou il est nécessaire de dilater les objets pieds des dunes
SE = strel([1 1;1 1]);
J = imdilate(DPP,SE);
DPP = J;

%Impression de l'image avec les numéros d'objet des pieds des dunes
DPP_Obj = 1:1:length(B_DPP);
TextDPP = string(DPP_Obj(:));

%calcul de l'angle pour le calcul pour jumeler les pieds des dunes à la
%ligne de crête
%Azimuth = Est la l'orientation de la ligne passant entre les points de
%début et fin des dunes. 
[Azimuth,Rumo] = Orientation_Dunes1(ID_Coord_All);
for i=1:length(Azimuth)
    %ID de la ligne de crête
    Angles(i,1) = ID_Coord_All(i,1);
    %Azimuth
    Angles(i,2) = Azimuth(i);
    %Nadir
    Angles(i,3) = Azimuth(i);
    %Angle perpendiculaire 1
    Angles(i,4) = Rumo(i)+90;    
    %Angle perpendiculaire 2
    Angles(i,5) = Rumo(i)-90;
end
%******************************************************************************************************************************************

%Jumeler les pieds des dunes avec les lignes de crête
%paramètres utilisées
%dist = distance de balayage des pixels (unité pixel)
dist = 1;
%pixel de début du pairage
debut = 1;
%Résolution - à chaque combien de pixels de la ligne de crête qu'il calcul
%le pied de dune jumeaux
res = 1;
%Maximum d'itérations dans le calcul
imax = 100;
%Ligne de pieds des dunes nord  - à noter nord et sud par rapport à matlab
[Obj_found_north,Obj_found_northX,Obj_found_northY, Vec_North_X, Vec_North_Y] = PiedsDunes1(B_DLC_Break,DLC_Break,Angles,L_DPP,DPP, imax, debut, res, dist);
%Ligne de pieds des dunes sud  - à noter nord et sud par rapport à matlab
[Obj_found_south,Obj_found_southX,Obj_found_southY, Vec_South_X, Vec_South_Y] = PiedsDunes2(B_DLC_Break,DLC_Break,Angles,L_DPP,DPP, imax, debut, res, dist);

%Fonction pour l'impression de l'objet de la ligne de crête avec ses pairs
%pieds des dunes - SEULEMENT les pixels touchés par le calcul précédent

%Selon le diagramme avec le modële opêrationnel - 2ème étape
TEST_Color_Dunes = zeros(size(L_DPP));
n1=1;
nf=length(B_DLC_Break);
m = TEST_Color_Dunes;
axis on;
axis image;
grid on;
DUNE_Provisoire = zeros(size(L_DPP));
SE = strel([1 1 1;1 1 1;1 1 1]);
for i=n1:nf
    X = Vec_South_X{i};
    Y = Vec_South_Y{i};
    [nt,mt] = size(Vec_South_X{i});
    for j=1:nt
        y1 = Vec_South_X{i}(j,1);
        x1 = Vec_South_Y{i}(j,1);
        y2 = Vec_South_X{i}(j,2);
        x2 = Vec_South_Y{i}(j,2);
        if ((y1==0)||(x1==0)||(y2==0)||(y2==0))
        else
        spacing =0.5;
        numSamples = ceil(sqrt((x2-x1)^2+(y2-y1)^2) / spacing);
        x = linspace(x1, x2, numSamples);
        y = linspace(y1, y2, numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0; sum(dxy, 2) == 0];
        finalxy = xy(~duplicateRows,:);
        finalx = finalxy(:, 1);
        finaly = finalxy(:, 2);
        for wt=1:length(finaly)
        DUNE_Provisoire(finaly(wt),finalx(wt)) = 1;
        end
        end
    end
     DUNE_Provisoire = imfill(DUNE_Provisoire,'holes');
     DUNE_Provisoire = imopen(DUNE_Provisoire,SE);
    [row col]=find(DUNE_Provisoire);
    for iw=1:length(row)
        if (TEST_Color_Dunes(row(iw),col(iw))~=0)
        else
           TEST_Color_Dunes(row(iw),col(iw)) = i; 
        end
    end
    DUNE_Provisoire = zeros(size(L_DPP));
    row = [];
    col = [];
end

hold on

DUNE_Provisoire = zeros(size(L_DPP));
for i=n1:nf
    X = Vec_North_X{i};
    Y = Vec_North_Y{i};
    obj = i;
    OBJ_Dune(i).DLC = i;
    [nt,mt] = size(Vec_North_X{i});
    for j=1:nt
        y1 = Vec_North_X{i}(j,1);
        x1 = Vec_North_Y{i}(j,1);
        y2 = Vec_North_X{i}(j,2);
        x2 = Vec_North_Y{i}(j,2);
        if ((y1==0)||(x1==0)||(y2==0)||(y2==0))
        else
        spacing =0.5;
        numSamples = ceil(sqrt((x2-x1)^2+(y2-y1)^2) / spacing);
        x = linspace(x1, x2, numSamples);
        y = linspace(y1, y2, numSamples);
        xy = round([x',y']);
        dxy = abs(diff(xy, 1));
        duplicateRows = [0; sum(dxy, 2) == 0];
        finalxy = xy(~duplicateRows,:);
        finalx = finalxy(:, 1);
        finaly = finalxy(:, 2);
        for wt=1:length(finaly)
        DUNE_Provisoire(finaly(wt),finalx(wt)) = 1;
        end
        end
    end
    DUNE_Provisoire = imfill(DUNE_Provisoire,'holes');
    DUNE_Provisoire = imopen(DUNE_Provisoire,SE);
    [row col]=find(DUNE_Provisoire);
    for iw=1:length(row)
        if (TEST_Color_Dunes(row(iw),col(iw))~=0)
            
        else
           TEST_Color_Dunes(row(iw),col(iw)) = i; 
        end
    end
    DUNE_Provisoire = zeros(size(L_DPP));
    row = [];
    col = [];
end

figure;imshow(TEST_Color_Dunes);

coordRefSysCode = 32187;
geotiffwrite('G15_19_1_100_OBJ_DUNES.tif',TEST_Color_Dunes, R1DTM,'CoordRefSysCode', coordRefSysCode)
figure;imshow(imfuse(DLC_Break,TEST_Color_Dunes)); 
hold on;plot(ID_Coord_All(:,2),ID_Coord_All(:,3),'r*');hold on
plot(ID_Coord_All(:,4),ID_Coord_All(:,5),'b*')
hold on;
Text = string(ID_Coord_All(:,1));
for i=1:length(ID_Coord_All(:,1))
text(ID_Coord_All(i,2),ID_Coord_All(i,3),Text(i),'Color','r','FontSize',14)
end
hold on
for i=1:length(ID_Coord_All(:,1))
text(ID_Coord_All(i,4),ID_Coord_All(i,5),Text(i),'Color','b','FontSize',14)
end
title('Lignes de crête avec les points de début et de fin'); hold on;
for i=2:2
    LigneCreteCoord = B_DLC_Break{i,:};
    [XLCC] = LigneCreteCoord(:,1);
    [YLCC] = LigneCreteCoord(:,2);
    plot(YLCC,XLCC,'r*');hold on;
        [XPD] = Obj_found_northX{i}(:);
        [YPD] = Obj_found_northY{i}(:);
        plot(YPD,XPD,'c*');hold on;
        [XPD] = Obj_found_southX{i}(:);
        [YPD] = Obj_found_southY{i}(:);
        plot(YPD,XPD,'y*');hold on;
    hold off

end

% --------------------------------------------------------------------------------------------------------------------------------------------------------
%METRIC DESCRIPTORS OF THE DUNE

%Enlever paires de coordonnées repétées dans les lignes de crête.
n_obj = length(B_DLC_Break);
for i=1:n_obj
    LigneCreteCoord = B_DLC_Break{i,:};
    A = unique(LigneCreteCoord,'rows');
    [Xo] = A(:,1);
    [Yo] = A(:,2);
    for j=1:length(A)
        Dune_Crest{i}(j,1) = Xo(j);
        Dune_Crest{i}(j,2) = Yo(j);
        Dune_Crest{i}(j,3) = DTM(Xo(j),Yo(j));
    end
end
%Length of the crest line (Lenght) and geodetic distance (Dc)
[Length,Dc] = Ditances(B_DLC_Break,ID_Coord_All);

%Calcul de la profondeur des dunes
n_obj = length(B_DLC_Break);
for i=1:n_obj
    Pc_moy{i} = mean(Dune_Crest{i}(:,3));
    Pc_med{i} = median(Dune_Crest{i}(:,3));
    Pc_max{i} = max(Dune_Crest{i}(:,3));
    Pc_min{i} = min(Dune_Crest{i}(:,3));
end

%Azimuth = Est la l'orientation de la ligne passant entre les points de
%début et fin des dunes. 
%AzimuthMigration = Direction d emigration des dunes. 
%calcul de l'orientation de chaq1ue objet dune;

AzimuthMigration = Azimuth - 90;
%Avoir une valeur d'Azimuth non négatif
for i=1:length(AzimuthMigration)
    if (AzimuthMigration(i)<0)
        AzimuthMigration(i) = 360+AzimuthMigration(i);
    else
    end
end
AzimuthMigration = AzimuthMigration';

%Calcul Width et hauteur des dunes
n_obj = length(B_DLC_Break);
for i=1:n_obj
    LigneCreteCoord = B_DLC_Break{i,:};
    [Xo] = LigneCreteCoord(:,1);
    [Yo] = LigneCreteCoord(:,2);
    if (isempty(Vec_North_X{i})==1 || isempty(Vec_North_Y{i})==1 || isempty(Vec_South_X{i})==1 || isempty(Vec_South_Y{i})==1 )
        Vec_North_X{i} = [1 1 1];
        Vec_North_Y{i} = [1 1 1];
        Vec_South_X{i} = [1 1 1];
        Vec_South_Y{i} = [1 1 1];
    else
    end
    for j=1:length(Xo)
        VXN = find(Vec_North_X{i}(:,3)==j);
        VYN = find(Vec_North_Y{i}(:,3)==j);
        VXS = find(Vec_South_X{i}(:,3)==j);
        VYS = find(Vec_South_Y{i}(:,3)==j);
        if (VXN~=0)&(VYN~=0)&(VXS~=0)&(VYS~=0)
            %Coordonnées ligne de crête
            Dune_Crest2{i}(j,1) = Xo(j);
            Dune_Crest2{i}(j,2) = Yo(j);
            Dune_Crest2{i}(j,3) = DTM(Xo(j),Yo(j));
            %Coordonnées pieds sud
            Dune_SouthTrough{i}(j,1) = Vec_South_X{i}(VXS,2);
            Dune_SouthTrough{i}(j,2) = Vec_South_Y{i}(VYS,2);
            Dune_SouthTrough{i}(j,3) = DTM(Vec_South_X{i}(VXS,2),Vec_South_Y{i}(VYS,2));
            %Coordonnées pieds nord
            Dune_NorthTrough{i}(j,1) = Vec_North_X{i}(VXN,2);
            Dune_NorthTrough{i}(j,2) = Vec_North_Y{i}(VYN,2);
            Dune_NorthTrough{i}(j,3) = DTM(Vec_North_X{i}(VXN,2),Vec_North_Y{i}(VYN,2));
            %Calcul des distances 2D - Width considéré dans l'article
            Width_Troughs{i}(j,1) = sqrt((Dune_SouthTrough{i}(j,1)-Dune_NorthTrough{i}(j,1))^2+(Dune_SouthTrough{i}(j,2)-Dune_NorthTrough{i}(j,2))^2);
            Width_CLST{i}(j,1) = sqrt((Dune_Crest2{i}(j,1)-Dune_SouthTrough{i}(j,1))^2+(Dune_Crest2{i}(j,2)-Dune_SouthTrough{i}(j,2))^2);
            Width_CLNT{i}(j,1)= sqrt((Dune_Crest2{i}(j,1)-Dune_NorthTrough{i}(j,1))^2+(Dune_Crest2{i}(j,2)-Dune_NorthTrough{i}(j,2))^2);
            %Calcul des distances 3D
            Distance_troughs{i}(j,1) = sqrt((Dune_SouthTrough{i}(j,1)-Dune_NorthTrough{i}(j,1))^2+(Dune_SouthTrough{i}(j,2)-Dune_NorthTrough{i}(j,2))^2+(Dune_SouthTrough{i}(j,3)-Dune_NorthTrough{i}(j,3))^2);
            Distance_CLST{i}(j,1) = sqrt((Dune_Crest2{i}(j,1)-Dune_SouthTrough{i}(j,1))^2+(Dune_Crest2{i}(j,2)-Dune_SouthTrough{i}(j,2))^2+(Dune_Crest2{i}(j,3)-Dune_SouthTrough{i}(j,3))^2);
            Distance_CLNT{i}(j,1) = sqrt((Dune_Crest2{i}(j,1)-Dune_NorthTrough{i}(j,1))^2+(Dune_Crest2{i}(j,2)-Dune_NorthTrough{i}(j,2))^2+(Dune_Crest2{i}(j,3)-Dune_NorthTrough{i}(j,3))^2);
            alpha{i}(j,1) = acosd((-(Distance_CLNT{i}(j,1))^2+(Distance_CLST{i}(j,1))^2+(Distance_troughs{i}(j,1))^2)/(2*(Distance_CLST{i}(j,1))*(Distance_troughs{i}(j,1))));
            H_alpha{i}(j,1) = sind(alpha{i}(j,1))*(Distance_CLST{i}(j,1));
            %%Calcul de des angles de pente (angle de stoss et lee side)
            alpha_South{i}(j,1) = acosd((sqrt((Dune_SouthTrough{i}(j,1)-Dune_Crest2{i}(j,1))^2+(Dune_SouthTrough{i}(j,2)-Dune_Crest2{i}(j,2))^2))/(Distance_CLST{i}(j,1)));
            alpha_North{i}(j,1) = acosd((sqrt((Dune_NorthTrough{i}(j,1)-Dune_Crest2{i}(j,1))^2+(Dune_NorthTrough{i}(j,2)-Dune_Crest2{i}(j,2))^2))/(Distance_CLNT{i}(j,1)));
        else
        
        end
    end
%     end
%     Hauteurs: Moyenne, Médianne, Maximale et minimale
            Hc_moy{i} = mean(H_alpha{i}(:,1));
            Hc_med{i} = median(H_alpha{i}(:,1));
            Hc_max{i} = max(H_alpha{i}(:,1));
            Hc_min{i} = min(H_alpha{i}(:,1));
%     Width Total: Moyen, Médian, Maximal et minimal
            Wd_moy{i} = mean(Width_Troughs{i}(:,1));
            Wd_med{i} = median(Width_Troughs{i}(:,1));
            Wd_max{i} = max(Width_Troughs{i}(:,1));
            Wd_min{i} = min(Width_Troughs{i}(:,1));
%     Width North: Moyen, Médian, Maximal et minimal
            WN_moy{i} = mean(Width_CLNT{i}(:,1));
            WN_med{i} = median(Width_CLNT{i}(:,1));
            WN_max{i} = max(Width_CLNT{i}(:,1));
            WN_min{i} = min(Width_CLNT{i}(:,1));
%     Width South: Moyen, Médian, Maximal et minimal
            WS_moy{i} = mean(Width_CLST{i}(:,1));
            WS_med{i} = median(Width_CLST{i}(:,1));
            WS_max{i} = max(Width_CLST{i}(:,1));
            WS_min{i} = min(Width_CLST{i}(:,1));
%     South angle: Moyen, Médian, Maximal et minimal
             alpha_South_moy{i} =  mean(alpha_South{i}(:,1));
             alpha_South_med{i} =  median(alpha_South{i}(:,1));
             alpha_South_max{i} =  max(alpha_South{i}(:,1));
             alpha_South_min{i} =  min(alpha_South{i}(:,1));
%     North angle: Moyen, Médian, Maximal et minimal
             alpha_North_moy{i} =  mean(alpha_North{i}(:,1));
             alpha_North_med{i} =  median(alpha_North{i}(:,1));
             alpha_North_max{i} =  max(alpha_North{i}(:,1));
             alpha_North_min{i} =  min(alpha_North{i}(:,1));
end

%Direction des migration des dunes
Om = AzimuthMigration;
%Crestline length
Lc = Length';
%Geodetic distance
Dc = Dc;
%Depth of the crest line
%Ordre - moy med max min
n_obj = length(B_DLC_Break);
for i=1:n_obj
    %Depth of the crest line
    Pc(i,1) = Pc_moy{i};
    Pc(i,2) = Pc_med{i};
    Pc(i,3) = Pc_max{i};
    Pc(i,4) = Pc_min{i};
    %Height of the crest line
    Hc(i,1) = Hc_moy{i};
    Hc(i,2) = Hc_med{i};
    Hc(i,3) = Hc_max{i};
    Hc(i,4) = Hc_min{i};
    
    %Width of the dune
    Wd(i,1) = Wd_moy{i};
    Wd(i,2) = Wd_med{i};
    Wd(i,3) = Wd_max{i};
    Wd(i,4) = Wd_min{i};
    
    %Width North
    WN(i,1) = WN_moy{i};
    WN(i,2) = WN_med{i};
    WN(i,3) = WN_max{i};
    WN(i,4) = WN_min{i};
    
    %Width South
    WS(i,1) = WS_moy{i};
    WS(i,2) = WS_med{i};
    WS(i,3) = WS_max{i};
    WS(i,4) = WS_min{i};
    
    %South ANGLE
    AS(i,1) = alpha_South_moy{i};
    AS(i,2) = alpha_South_med{i};
    AS(i,3) = alpha_South_max{i};
    AS(i,4) = alpha_South_min{i};

    %North ANGLE
    AN(i,1) = alpha_North_moy{i};
    AN(i,2) = alpha_North_med{i};
    AN(i,3) = alpha_North_max{i};
    AN(i,4) = alpha_North_min{i};
end

%Estimer un seul paramètre en considérant soit la mediane et en cas de
%nullité, utiliser la moyenne. 

n_obj = length(B_DLC_Break);
for i=1:n_obj
    %Depth of the crest line
    %Élimination de la moyenne et médianne en utilisant et conservant un seul paramètre
    if (Pc(i,2)==0)
        Pc1(i,1) = Pc(i,1); 
    else
        Pc1(i,1) = Pc(i,2); 
    end
    Pc1(i,2) = Pc(i,3);
    Pc1(i,3) = Pc(i,4);
    
    %Height of the crest line
    if (Hc(i,2)==0)
        Hc1(i,1) = Hc(i,1); 
    else
        Hc1(i,1) = Hc(i,2); 
    end
    Hc1(i,2) = Hc(i,3);
    Hc1(i,3) = Hc(i,4);
    
    %Width of the dune
    if (Wd(i,2)==0)
        Wd1(i,1) = Wd(i,1); 
    else
        Wd1(i,1) = Wd(i,2); 
    end
    Wd1(i,2) = Wd(i,3);
    Wd1(i,3) = Wd(i,4);
    
    %Width North
    if (WN(i,2)==0)
        WN1(i,1) = WN(i,1); 
    else
        WN1(i,1) = WN(i,2); 
    end
    WN1(i,2) = WN(i,3);
    WN1(i,3) = WN(i,4);

    %Width South
    if (WS(i,2)==0)
        WS1(i,1) = WS(i,1); 
    else
        WS1(i,1) = WS(i,2); 
    end
    WS1(i,2) = WS(i,3);
    WS1(i,3) = WS(i,4);
    
    %South ANGLE
    if (AS(i,2)==0)
        AS1(i,1) = AS(i,1); 
    else
        AS1(i,1) = AS(i,2); 
    end
    AS1(i,2) = AS(i,3);
    AS1(i,3) = AS(i,4);

    %North ANGLE
    if (AN(i,2)==0)
        AN1(i,1) = AN(i,1); 
    else
        AN1(i,1) = AN(i,2); 
    end
    AN1(i,2) = AN(i,3);
    AN1(i,3) = AN(i,4);
end

%Définir Quel côté est le lee et le stoss

n_obj = length(B_DLC_Break);
for i=1:n_obj
    if (AN1(i,1)>AS1(i,1))
        ALee(i,:) = AN1(i,:);
        WLee(i,:) = WN1(i,:);
        AStoss(i,:) = AS1(i,:);
        WStoss(i,:) = WS1(i,:);
    else

        ALee(i,:) = AS1(i,:);
        WLee(i,:) = WS1(i,:);
        AStoss(i,:) = AN1(i,:);
        WStoss(i,:) = WN1(i,:);
    end
end
%Index
% Sinousity Index
Si = Lc./Dc;
%Steepness Index
St = (Hc1(:,1))./(Wd1(:,1));
%Symmetry Index
Sy = (WStoss(:,1))./(WLee(:,1));
%Matrice avec tous les indices des dunes calculés. 
DunesParametres = [[1:1:n_obj]' Om Lc Dc Pc1(:,1:2) Hc1(:,1) Wd1(:,1:2) WStoss(:,2) WLee(:,2) AStoss(:,2) ALee(:,2) Si St Sy Hc1(:,2)];
[X_Duneincomplete,Y_Duneincomplete]=find(isnan(DunesParametres));
IndexX = unique(X_Duneincomplete);
[xtaille,ytaille] = size(IndexX);
for i=1:xtaille
    DunesParametres(IndexX(i),2:17)= [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
end
writematrix(DunesParametres,'G15_19_1_100_OBJ_DUNES.csv');

AFFICHE = figure;imshow(imfuse(DLC_Break,TEST_Color_Dunes,'Scaling','none'));
hold on;plot(ID_Coord_All(:,2),ID_Coord_All(:,3),'r*');hold on
plot(ID_Coord_All(:,4),ID_Coord_All(:,5),'b*')
hold on;
Text = string(ID_Coord_All(:,1));
for i=1:length(ID_Coord_All(:,1))
text(ID_Coord_All(i,2),ID_Coord_All(i,3),Text(i),'Color','r','FontSize',14)
end
hold on
for i=1:length(ID_Coord_All(:,1))
text(ID_Coord_All(i,4),ID_Coord_All(i,5),Text(i),'Color','b','FontSize',14)
end
title('Lignes de crête avec les objets dune'); hold on;
hold off
toc
