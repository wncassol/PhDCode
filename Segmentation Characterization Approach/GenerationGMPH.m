%Détection des Dunes
%Projet de Doctorat Willian Ney CASSOL
%Début de cette version du code 08 mai 2020
%Code adapté de la version développée par Vincent Dupont en geomorphon1

clear all
close all
format compact
format long g
clc
DTM_file = strcat(pwd,'/Export/G11/G11_310813_100.tif');

%Fonction Geomorphons créée sur Matlab par Vincent Dupont
%ParamÈtres des geomorphons par vecteur:
S=[10 10 10 10 10 10 10 10 15 20 10 15 15 15 30 20];
SKIP = [1 3 1 2 5 5 5 4 3 3 4 4 5 5 10 5];
F = [1 1 0.5 0.5 1 2 1.5 1.5 2 1 1 1 2 1.5 3 2];
D = [0 0 0 0 0 0 0 0 1 0.5 2 2 3 2 5 3];
[lg,kg]=size(S);

for i=1:kg
    filename = ['/resultat/GMPH/G11/G11_310813_100_gmph' num2str(i) '.tif'];
    gmph_file = strcat(pwd,filename);
    [gmph_out] = geomorphon1(DTM_file, S(i), SKIP(i), F(i), D(i), gmph_file);
end