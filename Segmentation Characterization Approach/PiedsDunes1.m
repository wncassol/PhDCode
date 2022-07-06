function [Obj_found_north,Obj_found_northX,Obj_found_northY, Vec_North_X, Vec_North_Y] = PiedsDunes1(B_DLC_Break,DLC_Break,Angles,L_DPP,DPP, imax, debut, res, dist)

n_obj = length(B_DLC_Break);

% Ligne de pieds des dunes nord - à noter nord et sud par rapport à matlab,
% laors inversé par apport aux vrai nord et sud
for i=1:n_obj
    %Coordonnées des pixels des objets lignes de crête
    C = B_DLC_Break{i,:};
    [Xo] = C(:,1);
    [Yo] = C(:,2);
    pos =1;
    for j=debut:res:length(Xo)
        Y0 = Yo(j);
        X0 = Xo(j);
            Xi = X0(1);
            Yi = Y0(1);
        iter=1;
        %itération jusqu'à trouver une ligne de pieds de dune jusqu'aux
        %paramètres prédéfinies
        while ((iter<=imax))
            Yp = Y0 + dist*sind(Angles(i,5));
            Xp = X0 + dist*cosd(Angles(i,5));
            iter = iter+1;
            Xpr = round(Xp);
            Ypr = round(Yp);
            [k,l] = size(L_DPP);
            if ((Xpr>k) || (Ypr>l) || (Xpr==0) || (Ypr==0))
                break;
            else
                if (DPP(Xpr,Ypr)>0)
                %Enregistrement de l'objets pied de dune trouvé pour cette
                %ligne de crête
                Obj_found_north{i}(pos,1) = L_DPP(Xpr,Ypr);
                %enregistrement des coordonnées X et Y touchées par le
                %calcul de jumeleage
                Obj_found_northX{i}(pos,1) = Xpr;
                Obj_found_northY{i}(pos,1) = Ypr;
                Vec_North_X{i}(pos,1:3) = [Xi Xpr j];
                Vec_North_Y{i}(pos,1:3) = [Yi Ypr j];
                pos=pos+1;
                break;
                elseif (DLC_Break(Xpr,Ypr)>0)
                Obj_found_north{i}(pos,1) = 0;
                else
                X0 = Xpr;
                Y0 = Ypr;
                end
            end
        end
   end
end
 

if (length(Obj_found_north)<n_obj)
    Obj_found_north{n_obj}(1,1) = 0;
    Vec_North_X{n_obj}(1,1:2) = [0 0];
    Vec_North_Y{n_obj}(1,1:2) = [0 0];
end
