function [gmph_out] = geomorphon1(mnt_file, search, skip, flat, dist, gmph_file)
% Prototype de fonction: [nuage_class, gmph_out] = geomorphon(nuage_file, mnt_file, search, skip, flat, dist, nuage_class_file, gmph_file)
%
% R�le de la fonction : 
%   Cette fonction permet de classifier les 10 classes g�omorphologiques 
%   d'un MNT (Mod�le Num�rique de Terrain). Pour se faire, l'outil 
%   r.geomorphon disponible dans la librairie GRASS GIS est utilis�. Le 
%   premier r�sultat de cette classification est une grille matricielle 
%   appel�e 'g�omorphon'. Le deuxi�me r�sultat est une nuage de points 
%   classifi� de la valeur du g�omorphon. Autrement dit, la valeur de la 
%   grille du g�omorphon est associ�e aux points inclus dans la cellule 
%   classifi�e de la grille matricielle des g�omorphons. 
%   De plus, ajoutons que le logiciel SAGA GIS est utilis� pour pouvoir 
%   lire le nuage de points et associer les valeurs de la grille des 
%   g�omorphons au nuage de points.
%
%   Pr�conditions : 
%   1) Pour utiliser cette fonction. L'utilisateur doit d�finir
%   manuellement dans ce code les chemins d'acc�s en absolu des fichiers 
%   saga_cmd.exe et grass78.bat install� sur son ordinateur pour que la 
%   fonction Matlab puisse utiliser SAGA GIS et GRASS GIS. Pour se faire, 
%   modifiez les variables SAGA_cmd et GRASSBAT.
%   
%   2) Le syst�me de coordonn�es du nuage de points et du MNT doit �tre le
%      m�me.
%
%   Pour plus d'information sur la fonction r.geomorphon et les logiciels
%   GRASS GIS et SAGA GIS, r�f�rez vous aux liens suivants : 
%
%   GRASS GIS : https://grass.osgeo.org/
%   GRASS GIS fonction r.geomorphon : 
%   https://grass.osgeo.org/grass78/manuals/r.geomorphon.html
%
%   SAGA GIS : http://www.saga-gis.org/saga_tool_doc/6.3.0/index.html
%   
% Arguments d'entr�e : 
%
%  nuage_file : Chemin d'acc�s absolu du fichier xyz du nuage de points. 
%               Le fichier doit avoir 3 colonnes sans en-t�te s�par�es par 
%               un espace. (type character)
%               Exemple : 'C:/Projet/data/nuage.xyz'
%
%  mnt_file   : Chemin d'acc�s absolu du fichier GeoTiff du MNT cr�� � 
%               partir du nuage de points 'nuage_file'. Le syst�me de 
%               r�f�rence cartographique doit �tre d�fini et �tre le m�me 
%               que celui des coordonn�es du nuage de points. (type
%               character)
%               Exemple : 'C:/Projet/data/mnt.tif'
%
%  search, skip, flat, dist : Param�tres pour appliquer la fonction 
%                             r.geomorphon(type entier). 
%                             Voir la documentation : 
%                 https://grass.osgeo.org/grass78/manuals/r.geomorphon.html
%
%
% Param�tres de sortie:
%
% nuage_class_file : Chemin d'acc�s absolu du fichier txt du nuage de 
%                    points classifi�. Le fichier de sortie a 3 colonnes 
%                    sans en-t�te s�par�es par un espace.(type character)
%                    Exemple : 'C:/Projet/resultat/nuage_class.txt'
%
% gmph_file   : Chemin d'acc�s absolu du fichier GeoTiff du g�omorphon 
%               Le syst�me de r�f�rence cartographique est le m�me que 
%               celui du MNT.(type character)
%               Exemple : 'C:/Projet/resultat/geomorph.tif'
%
% nuage_class : Matrice Matlab Nx4 contenant les coordonn�es des points et
%               la valeur de la classe de g�omorphon associ�e. Ces valeurs
%               vont de 1 � 10. S'il n'y a aucune classe, la valeur -99999
%               est inscrite.
%
% gmph_out    : Structure contenant des �l�ments Matlab produits par les
%               fonctions geotiffread() et geotiffinfo().
%               gmph_out.value : Grille matricielle des valeurs du
%                                g�omorphon. S'il n'y pas de classe, la
%                                valeur NaN est inscrite.
%               gmph_out.refmat : Voir la documentation de 
%                                 geotiffread() et la variable refmat.
% 
%               gmph_out.info   : Structure des propri�t�s de le la grille
%                                 GeoTIFF file. Voir la documentation de la
%                                 fonction geotiffinfo().
%
% Exemple pour afficher la grille de g�omorphon dans Matlab : 
% figure;
% mapshow(gmph_out.value,gmph_out.info.SpatialRef,'DisplayType','surface')
% axis image off;
%
%
% Auteur : Vincent Dupont
% D�partement des sciences g�omatiques - Universit� Laval
% Date: 8 mai 2020
%
% -----------------------------------------------------------------------

% Chemin d'acc�s � SAGA GIS
SAGA_cmd = strcat('"',"C:/Program Files/SAGA/saga_cmd.exe",'"');

% Chemin d'acc�s vers le fichier grass78.bat
GRASSBAT = strcat('"',"C:/Program Files/GRASS GIS 7.8/grass78.bat",'"');

% Variables pour les tests (voir la fonction run() � la fin du code)
affiche = false;
act = true;

% Configurer les chemins d'acc�s
% nuage_file = replace(nuage_file,"\","/"); %input data
mnt_file = replace(mnt_file,"\","/"); %input data
% nuage_class_file = replace(nuage_class_file,"\","/"); %output data
gmph_file = replace(gmph_file,"\","/"); %output data

% Initialiser certaines variables utilitaires
% [nuage.filepath,nuage.name,nuage.ext] = fileparts(nuage_file);

[mnt.filepath,mnt.name,mnt.ext] = fileparts(mnt_file);
mnt.info = geotiffinfo(mnt_file);
mnt.epsg = num2str(mnt.info.GeoTIFFCodes.PCS);

[gmph.filepath,gmph.name,gmph.ext] = fileparts(gmph_file);


% �tape 1 : Initialiser GRASS GIS
% Cr�er le fichier GRASSDB pour GRASS GIS
GRASSDB = gmph.filepath+"/grassdb/";
LOCATION = "mylocation/";
MAPSET = "PERMANENT";
GRASSEXEC = GRASSBAT + " " + GRASSDB+LOCATION+MAPSET + " --exec";

cmd = [GRASSBAT, "-e", "-c", "EPSG:"+mnt.epsg, GRASSDB+LOCATION];

% Cr�er le dossier des GRASSDB s'il n'existe pas
if ~exist(GRASSDB, "dir")
    mkdir(convertStringsToChars(GRASSDB));
end

% Cr�er la base de donn�es GRASS si elle n'existe pas
if ~exist(GRASSDB+LOCATION, "dir")
    status = run(cmd, affiche, act);
else
    rmdir(convertStringsToChars(GRASSDB),'s');
    status = run(cmd, affiche, act);
end


% �tape 2 : Appliquer le geomorphon avec la fonction r.geomorphon
% Importer la grille du MNT dans GRASS
cmd = [GRASSEXEC, "r.external",...
       strcat("input=",mnt_file),...
       "output=" + mnt.name,...
       "--overwrite"];

status = run(cmd, affiche, act);

% Definir la zone de calcul
cmd = [GRASSEXEC, "g.region",...
       "raster=" + mnt.name,...
       "--overwrite"];

status = run(cmd, affiche, act);

% Appliquer le geomorphon sur le MNT
cmd = [GRASSEXEC, "r.geomorphon",...
       "elevation=" + mnt.name,...
       "search=" + num2str(search),...
       "skip=" + num2str(skip),...
       "flat=" + num2str(flat),...
       "dist=" + num2str(dist),...
       "forms=" + mnt.name + "_geomorph",...
       "--overwrite"];

status = run(cmd, affiche, act);

% Exporter la grille du g�omorphon en format GeoTiff
cmd = [GRASSEXEC, "r.out.gdal",...
       "input=" + mnt.name + "_geomorph",...
       strcat("output=",gmph_file),...
       "--overwrite"];

status = run(cmd, affiche, act);

% Exporter le geomorphon en format grille matricielle matlab
[gmph_out.value, gmph_out.refmat] = geotiffread(gmph_file);
gmph_out.info = geotiffinfo(gmph_file);

gmph_out.value = single(gmph_out.value);
gmph_out.value(gmph_out.value==255) = NaN;


% �tape 3 : Classifier le nuage de points avec les classes du g�omorphon 
% Lire le fichier xyz du nuage de points
% SAGA : Import Export -> Shapes -> Import Point Cloud from Text File
% % % cmd = [SAGA_cmd, "io_shapes", '16',...
% % %        strcat('-FILE=',nuage_file),...
% % %        "-SEPARATOR=1",...
% % %        "-SKIP_HEADER=0",...
% % %        "-XFIELD=1",...
% % %        "-YFIELD=2",...
% % %        "-ZFIELD=3",...
% % %        "-POINTS="+GRASSDB+nuage.name+".sg-pts-z"];
% % % status = run(cmd, affiche, act);
% % % 
% % % 
% % % % Assigner les attributs de la grille du g�omorphon aux points du nuage
% % % % SAGA : Shapes -> Grid Tools -> Add Grid Values to Shapes
% % % cmd = [SAGA_cmd, "shapes_grid", "1",...
% % %        "-SHAPES=" + GRASSDB+nuage.name+".sg-pts-z",...
% % %        "-GRIDS=" + gmph_file,...
% % %        "-RESAMPLING=0"];
% % % status = run(cmd, affiche, act);
% % % 
% % % % Exporter en format txt le nuage de points classifi�
% % % % SAGA Import Export --> Shapes --> Export Point Cloud to Text File
% % % cmd = [SAGA_cmd, "io_shapes", "18",...
% % %        "-POINTS=" + GRASSDB+nuage.name+".sg-pts-z",...
% % %        "-FILE=" + nuage_class_file,...
% % %        "-WRITE_HEADER=0",...
% % %        "-FIELDSEP=1",...
% % %        "-FIELDS=" + strcat('"',"1;2;3;4",'"'),...
% % %        "-PRECISIONS=" + strcat('"',"3;3;3;-1",'"')];
% % % status = run(cmd, affiche, act);
% % % 
% % % % Lire le fichier du nuage de points classifi�
% % % nuage_class = load(nuage_class_file);

% �tape 4 : D�truire le dossier GRASSDB contenant des fichiers temporaires
if exist(GRASSDB, "dir")
   rmdir(convertStringsToChars(GRASSDB),'s');
end


% Fonction pour lancer les commandes et faire des tests.
% affiche : Permet d'afficher dans la console les traitements SAGA et GRASS
% act : Ex�cuter les commandes SAGA et GRASS (true) ou affiche le string de
% la commande dans la console sans ex�cuter les commandes.(Si affiche true) 
function status = run(cmd, affiche, act)

status = 0;
if(act)
    if(affiche)
        disp(join(cmd));
        status = system(join(cmd));
    else
        [status,~] = system(join(cmd));
    end
else
    disp(join(cmd));
end

end % fin de run()

end % fin de geomorphon()