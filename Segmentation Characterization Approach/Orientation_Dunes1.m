function [Azimute,Rumo] = Orientation_Dunes(ID_Coord_All)

%Projet de Doctorat Willian Ney CASSOL
%Début de cette version du code 01 octobre 2021
%Cette fonction calcule les gisements, azimuts ou l'orientation des dunes
%considérant leur début et fin. 
%ID_Coord_All = ID Ligne Crete, Xd, Yd, Xf, Yf
ID = ID_Coord_All(:,1);
Xd = ID_Coord_All(:,2);
Yd = ID_Coord_All(:,3);
Xf = ID_Coord_All(:,4);
Yf = ID_Coord_All(:,5);
%Geodetic distance
% Dc = sqrt((Xf-Xd).^2+(Yf-Yd).^2);

for i=1:length(ID)
    
    DeltaY(i) = Yf(i)-Yd(i);
    DeltaX(i) = Xf(i)-Xd(i);
    Rumo(i) = atan((DeltaX(i))/(DeltaY(i)));
    Rumo(i)=rad2deg(Rumo(i));
    
    if (DeltaX(i)>0) && (DeltaY(i)>0)
        Azimute(i) = Rumo(i);        
        
    elseif (DeltaX(i)>0) && (DeltaY(i)<0)
        Azimute(i) = 180-Rumo(i); 
        
    elseif (DeltaX(i)<0) && (DeltaY(i)<0)
        Azimute(i) = 180+Rumo(i); 
        
    else
        Azimute(i) = 360-Rumo(i); 
        
    end
    

    %puisque l'image est orientée vers le bas, il faut l'orienter vers le
    %Nord, le haut dans matlab
    
    if (Azimute(i)>0) && (Azimute(i)<180)
        Azimute(i) = 180-Azimute(i); 
    elseif (Azimute(i)>=180) && (Azimute(i)<=360)
        Azimute(i) = Azimute(i)-180; 
    end
    

end
