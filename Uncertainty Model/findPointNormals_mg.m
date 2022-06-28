function [normales, courbure] = findPointNormals_mg(points, voisins, p)
% Prototype de fonction: [normales, courbure] = findPointNormals_mg(points, voisins, p)
%
% Rôle de la fonction : 
%
% Cette fonction permet de calculer le vecteur normal et la courbure locale 
% pour chacun des points du nuage en utilisant un ACP robuste. Cette
% dernière exploite le concepte de la moyenne généralisée (mg).
% 
% Arguments en entrée : 
%
%  points : Matrice nx3 ou chaque ligne représente un point (x,y,z).
%  voisins : Le nombre de voisins considérés pour calculer la normale et la
%            courbure (typiquement 20 voisins)
%  p : Paramètre p de la moyenne généralisée (typiquement 0.1)
%
%
% Paramètres de sortie:
%
%   normales : matrice nx3 contenant les normales locales unitaires (nx,ny,nz)
%   courbure :  vecteur nx1 contenant les courbures locales. 
%
% Auteur : Vincent Dupont
%
% NOTE IMPORTANTE : Ce code est une combinaison du code [1] en ajoutant les 
%                   éléments de ACP robuste que l'on retrouve dans [2].
%
% Date: 9 mai 2019
%
% Références
% [1]
%   Zachary Taylor
%   zacharyjeremytaylor@gmail.com
%   http://www.zjtaylor.com
%
%   http://pointclouds.org/documentation/tutorials/normal_estimation.php
%   This code was used in generating the results for the journal paper
%   Multi-modal sensor calibration using a gradient orientation measure 
%   http://www.zjtaylor.com/welcome/download_pdf?pdf=JFR2013.pdf
%   Code initial : https://www.mathworks.com/matlabcentral/fileexchange/48111-find-3d-normals-and-curvature 
%   (consulté le 9 mai 2019)
%
% [2] OH et N. Kwak, « Robust PCA », in Advances in Principal Component Analysis, 2017, 
%     chap. 4, p. 71-98, isbn : 9789811067037. doi : 10.1007/978-981-10-6704-4
%
% -----------------------------------------------------------------------

% Nombre de points
nb_pts = size(points,1);


% Objet de recherche Kdtree
kdtreeobj = KDTreeSearcher(points,'distance','euclidean');

% Extraction des indices des k plus proche voisins pour tous les points
n = knnsearch(kdtreeobj, points,'k',(voisins+1));

% Calcul de la moyenne généralisée pour chaque ensemble point et ses k voisins
mg_mat = zeros(nb_pts, 3);

for i = 1:nb_pts
    mg_mat(i,:) = generalized_mean(points(n(i,:)',:), p);  
end

% Élimination de la colonne 1 représentant l'indice du point où on
% recherche ses voisins
n = n(:,2:end);


% Calcul du delta X,Y,Z entre la moyenne généralisée du groupe des k voisins
% et tous les voisins
p0 = repmat(mg_mat(:,1:3),voisins,1) - points(n(:),1:3);

% Réorganisation de p pour pouvoir calculer la matrice de variance-covariance
% totale. Ici on reshape les dimensions sur la 3e dimension d'une matrice.
p = reshape(p0, size(points,1), voisins,3);

% Calcul de la matrice de variance-convariances générale. Chaque ligne représente les
% 6 éléments de la matrice de variance-covariance d'un point par rapport à ses k voisins 
C = zeros(size(points,1),6);
C(:,1) = sum(p(:,:,1).*p(:,:,1),2);
C(:,2) = sum(p(:,:,1).*p(:,:,2),2);
C(:,3) = sum(p(:,:,1).*p(:,:,3),2);
C(:,4) = sum(p(:,:,2).*p(:,:,2),2);
C(:,5) = sum(p(:,:,2).*p(:,:,3),2);
C(:,6) = sum(p(:,:,3).*p(:,:,3),2);
C = C ./ voisins;

% Calcul de la normale et de la courbure locale pour chaque point
normales = zeros(size(points));
courbure = zeros(size(points,1),1);

for i = 1:(size(points,1))
    
    % Sous matrice de variance-covariance des k voisins
    Cmat = [C(i,1) C(i,2) C(i,3);...
        C(i,2) C(i,4) C(i,5);...
        C(i,3) C(i,5) C(i,6)];
    
    % Calcul des vecteurs propres (v) et des valeurs propres (d) de la matrice Cmat
    [v,d] = eig(Cmat);
    
    % Extraction de la plus petite valeur propre
    d = diag(d);
    [lambda,k] = min(d);
    
    % Vecteur propre associé à la plus petite valeur propre.
	% Il correspond au vecteur normal local 
    vn = v(:,k)';
    
	% Si la 3e composante est négative, on inverse le vecteur pour que celui-ci pointe vers le haut
    if(vn(3) < 0)
        vn = vn.*-1;
    end
    
    % Enregistre la normale
    normales(i,:) = vn;
    
    % Enregistre la courbure
    courbure(i) = lambda / sum(d);
end



function [mg] = generalized_mean(pts, p)
% Prototype de fonction: [mg] = generalized_mean(pts, p)
%
% Rôle de la fonction : 
%
% Cette fonction permet de calculer la moyenne généralisée pour une nuage de
% points x,y,z. 
% 
% Arguments d'entrée : 
%
%  pts : Matrice nx3 ou chaque ligne représente un point.
%  p : Paramètre p de la moyenne généralisée (typiquement 0.1)
%
% Paramètres de sortie:
%
%  mg : Vecteur 1x3 représentant la moyenne généralisée des coordonnées x,y,z
%       obtenue après itération.
%
% Code réécrit par : Vincent Dupont
%
% Référence, algorithme et code original : 
% OH et N. Kwak, « Robust PCA », in Advances in Principal Component Analysis, 2017, 
% chap. 4, p. 71-98, isbn : 9789811067037. doi : 10.1007/978-981-10-6704-4
%
% Date: 9 mai 2019
%
% -----------------------------------------------------------------------

% Nombre de points
nb_pts = size(pts,1); %[nx3]

% Calcul de la moyenne arithmétique
moy = mean(pts,1); %[1x3]
moy_mat = repmat(moy, nb_pts,1); %[nx3] copier n fois moy

% Soustraire la moyenne arithmétique des coordonnées des points
delta_ptsMean = pts - moy_mat; % [nx3] - [nx3] = [nx3]

% Calcul de la norme de chaque vecteur d'écart (pts-moy) 
% Ensuite, on met à la puissance p chaque norme et on fait la somme
objMG = sum(diag(delta_ptsMean * delta_ptsMean').^p); % [1x1]

% Variable temporaire pour les moyennes
mg_i = moy;

% Critères d'arrêt du processus itératif
maxItr = 50;
seuil = 0.01;
flag = 0;

j = 1;
while (flag == 0)
    
    % variable pour le critère d'arrêt
    objMG_avant = objMG;
    mg_i_avant = mg_i;
    
    % Calcul des écarts entre le point et moyenne mg_i
    mg_mat = repmat(mg_i, nb_pts,1);
    delta_ptsMg = pts - mg_mat; %[nx3]
    
    % Calcul de alphas pour chaque composante
    alphas = diag(delta_ptsMg * delta_ptsMg').^(p-1); %[nx1]
    
    % Ici les 3 alphas multiplient chacune des coordonnées
    % [pts'*alphas] fait la somme pondérée
    mg_i = ( (pts'*alphas) / sum(alphas) )'; %[1x3] 
    
    % Avec le nouveau mg_i, tu peux recalculer l'écart entre les points et ce mg_i
    mg_mat = repmat(mg_i, nb_pts,1);
    delta_ptsMg = pts - mg_mat; %[nx3]
    
    % Calcul de la fonction à optimiser avec le new mg_i
    objMG = sum(diag(delta_ptsMg * delta_ptsMg').^p);
    
    % Écart pour le critère d'arrêt
    diff_mg = mg_i - mg_i_avant;
    diff_mg_norm = sqrt(diff_mg*diff_mg');
    
    % Critère d'arrêt
    if (j >= maxItr)
        flag = 1;
    elseif diff_mg_norm/sqrt(mg_i*mg_i')*100 < seuil
        flag = 1;
    elseif objMG >= objMG_avant
        flag = 1;
        mg_i = mg_i_avant;
    else
        j = j + 1;
    end
end

% Retour de la moyenne généralisée
mg = mg_i;

end




end