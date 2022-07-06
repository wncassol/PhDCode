function [Obj_found_south,Obj_found_southX,Obj_found_southY, Vec_South_X, Vec_South_Y] = PiedsDunes2(B_DLC_Break,DLC_Break,Angles,L_DPP,DPP, imax, debut, res, dist)

n_obj = length(B_DLC_Break);
% Ligne de pieds des dunes sud - à noter nord et sud par rapport à matlab,
% laors inversé par apport aux vrai nord et sud
Azimuth_Median2 = median(Angles(:,4));

for i=1:n_obj
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
        while ((iter<=imax))
            Yp = Y0 + dist*sind(Angles(i,4));
            Xp = X0 + dist*cosd(Angles(i,4));
            Xpr = round(Xp);
            Ypr = round(Yp);
            [k,l] = size(L_DPP);
            if ((Xpr>k) || (Ypr>l) || (Xpr==0) || (Ypr==0))
                break;
            else
                if (DPP(Xpr,Ypr)~=0)
                Obj_found_south{i}(pos,1) = L_DPP(Xpr,Ypr);
                Obj_found_southX{i}(pos,1) = Xpr;
                Obj_found_southY{i}(pos,1) = Ypr;
                Vec_South_X{i}(pos,1:3) = [Xi Xpr j];
                Vec_South_Y{i}(pos,1:3) = [Yi Ypr j];
                pos=pos+1;
                break;
                elseif (DLC_Break(Xpr,Ypr)>0)
                Obj_found_south{i}(pos,1) = 0;
                pos=pos+1;
                else
                X0 = Xpr;
                Y0 = Ypr;
                end
            end
        iter = iter+1;
        end
    end
end

if (length(Obj_found_south)<n_obj)
    Obj_found_south{n_obj}(1,1) = 0;
    Vec_South_X{n_obj}(1,1:2) = [0 0];
    Vec_South_Y{n_obj}(1,1:2) = [0 0];
end
