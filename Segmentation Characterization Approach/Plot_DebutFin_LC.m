function [Coord] = Plot_DebutFin_LC(B_DLC_Break,DPP,DLC_Break,area)
%Cette Fonction d�velopp�e dans le cadre du doctorat de Willian Ney Cassol
%plot le d�but et fin des objets lignes de cr�te identifi�s automatiquement
%par matlab. - Fonction peut-�tre obselete


n_obj = length(B_DLC_Break);
% n_obj = 20;
WT = imfuse(DLC_Break,DPP);

figure;imshow(WT);hold on;

for i=1:n_obj
    m = area(i);
    Vector = B_DLC_Break{i};
    [X_i] = Vector(1,1);
    [Y_i] = Vector(1,2);
    [X_f] = Vector(m,1);
    [Y_f] = Vector(m,2);
    WT = imfuse(DLC_Break,DPP);
    plot(Y_i,X_i,'r*');plot(Y_f,X_f,'b*');hold on
end
Coord = [X_i, Y_i, X_f, Y_f];
