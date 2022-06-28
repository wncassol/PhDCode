clear 
clc
close all
%Ce code répresente le chapitre 3 du mémoire.
%Legende
% s sigma - écart-type ou covariance
% S Sigma - Variance (Matrice)
% J - Derivées partielles
% d delta 
d2r = pi/180;
%**************************************************************************
% Simulateur des dunes - voir script simulateur dunes
% % Profil Fond plat
% Dune_Profil = load('Dune1.txt');
% Wp = load('rho1.txt');

% Profil au long du current sea-bottom
Dune_Profil = load('Dune1.txt');
Wp = load('rho1.txt');

% Profil Transversal
% Dune_Profil = load('Dune3.txt');
% Wp = load('rho3.txt');


rho = Wp(1,:);
gama = Wp(2,:);
SVP_SIMULE = 1471.4;
X = 0;
Wpente = pente_app(Dune_Profil);
%Entrée des données mésurées-----------------------------------------------
%Entrés(unités: m, s, m/s, rad, rad/s)
%Écarts-type des coordonnées
sdxn = 0.05;
sdyn = 0.05;
sdzn = 0.08;
%Latence entre LiDAR et recépteur GNSS
dTlp = 0.000;
%Latence entre Inertielle et LiDAR
dTli = 0.000;
%Paramètres de montage
%Lever-arms
ax = 0;
ay = 0;
az = 0;
abi = [ax;ay;az];
%Écarts-type des bras de levier
sdax = 0.015;
sday = 0.015;
sdaz = 0.015;
%Angles de Boresight
psib = 0*d2r;
thetab = 0*d2r;
phib = 0*d2r;
%Écarts-type des angles de boresight
sdphib = 0.05*d2r;
sdthetab = 0.05*d2r;
sdpsib = 0.05*d2r;
%Écarts-types des données fournies par le MBES
sdrho = 0.0;
sdgama =  0.002*d2r;
stdMBES = 0.000020;
scel = 0.025;
sCEL = [scel;sdgama];

[l,k] = size(gama);
for i=1:k
if any(isnan(Wpente(i)))
    Wpente(i)= Wpente(i-1);
else
end
alpha(i) = Wpente(i)*d2r;
%Paramètres du plan associés à l'inclinaison des dunes (alpha)
S = [0 sin(alpha(i)) cos(alpha(i))];
%**************************************************************************
xn = 0;
yn = 0;
zn = 0;
Pn = [xn;yn;zn];
%Vitesses Linéaires de la plateforme
VN = 0;
VE = 0;
VD = 0;
%Centrale inertielle (phi, theta, psi, w1, w2, w3)
%Angles d'attitude
phi = 0*d2r;
theta = 0*d2r;
psi = 0*d2r;
%Écarts-type des angles d'attitude
sdphi = 0.01*d2r;
sdtheta = 0.01*d2r;
sdpsi = 0.03*d2r;
%Vitesses Angulaires
w1 = 0*d2r;
w2 = 0*d2r;
w3 = 0*d2r;
Omega = [0 -w3 w2;w3 0 -w1;-w2 w1 0];
%--------------------------------------------------------------------------
%Matrice Variance-Covariance pour la position(Pn) -->Pn=Pn(Xn,Yn,Zn,t)
CoVar_xi = [sdxn^2 0 0 0;
         0 sdyn^2 0 0;
         0 0 sdzn^2 0;
         0 0 0 dTlp^2];
%Derivées partielles de Pn (Matrice des derivées Position)
JPn = [1 0 0 VN;
       0 1 0 VE;
       0 0 1 VD];
%Matrice Variance-Covariance Positionnement
CoVar_Pn = JPn*CoVar_xi*JPn';
%--------------------------------------------------------------------------
%Matrice Variance-Covariance lever-arms(an) -->zeta=(phi,theta,psi,t,ax,ay,az)
%Derivées partielles de an (Matrice des derivées partielles Lever-arms)
Jan = [R1(psi)*R2(theta)*JR3(phi)*abi R1(psi)*JR2(theta)*R3(phi)*abi JR1(psi)*R2(theta)*R3(phi)*abi Cbn(phi,theta,psi)*Omega*abi Cbn(phi,theta,psi)];
%MVC des paramètres de l'effet des lever-arms
CoVar_zeta = [sdphi^2 0 0 0 0 0 0;
       0 sdtheta^2 0 0 0 0 0;
       0 0 sdpsi^2 0 0 0 0;
       0 0 0 dTli^2 0 0 0;
       0 0 0 0 sdax^2 0 0;
       0 0 0 0 0 sday^2 0;
       0 0 0 0 0 0 sdaz^2];
CoVar_an = Jan*CoVar_zeta*Jan';
ubs = [0;sind(gama(i));cosd(gama(i))];
rbs = rho(i)*ubs;
Jrbs_rho =[0; sind(gama(i)); cosd(gama(i))];
Jrbs_gama =rho(i)*[0;cosd(gama(i));-sind(gama(i))];
%Écarts-types des données fournies par le MBES
wtt = (rho(i)/SVP_SIMULE)*2;
sdrho = [1/2*SVP_SIMULE 1/2*wtt]*[stdMBES^2 0;0 sCEL(1)^2]*[1/2*SVP_SIMULE; 1/2*wtt];
Jrnphi = R1(psi)*R2(theta)*JR3(phi)*Cbn(phib,thetab,psib)*rbs;
Jrntheta = R1(psi)*JR2(theta)*R3(phi)*Cbn(phib,thetab,psib)*rbs;
Jrnpsi = JR1(psi)*R2(theta)*R3(phi)*Cbn(phib,thetab,psib)*rbs;
Jrnt = Cbn(phi,theta,psi)*Omega*Cbn(phib,thetab,psib)*rbs;
Jrnphib = Cbn(phi,theta,psi)*R1(psib)*R2(thetab)*JR3(phib)*rbs;
Jrnthetab = Cbn(phi,theta,psi)*R1(psib)*JR2(thetab)*R3(phib)*rbs;
Jrnpsib = Cbn(phi,theta,psi)*JR1(psib)*R2(thetab)*R3(phib)*rbs;
Jrnrho = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*Jrbs_rho;
Jrngama = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*Jrbs_gama;
Jrn = [Jrnphi Jrntheta Jrnpsi Jrnt Jrnphib Jrnthetab Jrnpsib Jrnrho Jrngama];
%Dérivées partielles de Xn 
JXnphi = Jrn(:,1) + Jan(:,1); 
JXntheta = Jrn(:,2) + Jan(:,2); 
JXnpsi = Jrn(:,3) + Jan(:,3); 
JXnt = Jrn(:,4) + Jan(:,4); 
%Calcul de la matrice variance covariance pour rn
CoVar_chi = [sdphi^2 0 0 0 0 0 0 0 0;
       0 sdtheta^2 0 0 0 0 0 0 0;
       0 0 sdpsi^2 0 0 0 0 0 0;
       0 0 0 dTli^2 0 0 0 0 0;
       0 0 0 0 sdphib^2 0 0 0 0;
       0 0 0 0 0 sdthetab^2 0 0 0;
       0 0 0 0 0 0 sdpsib^2 0 0;
       0 0 0 0 0 0 0 sdrho^2 0;
       0 0 0 0 0 0 0 0 sdgama^2];
CoVar_rn = Jrn*CoVar_chi*Jrn';
den1 = S*Jrn(:,8);
JXn_chi = [JXnphi JXntheta JXnpsi JXnt Jrn(:,5) Jrn(:,6) Jrn(:,7) Jrn(:,8) Jrn(:,9)];
JXn_chi1 = [JXnphi JXntheta JXnpsi JXnt Jrn(:,5) Jrn(:,6) Jrn(:,7) Jrn(:,9)];
dchi1 = [sdphi; sdtheta; sdpsi; dTli; sdphib; sdthetab; sdpsib; sdgama];
num1 = S*JXn_chi1*dchi1;
drho_i =-(num1/den1);
Jsnphi = R1(psi)*R2(theta)*JR3(phi)*Cbn(phib,thetab,psib)*ubs;
Jsntheta = R1(psi)*JR2(theta)*R3(phi)*Cbn(phib,thetab,psib)*ubs;
Jsnpsi = JR1(psi)*R2(theta)*R3(phi)*Cbn(phib,thetab,psib)*ubs;
Jsnt = Cbn(phi,theta,psi)*Omega*Cbn(phib,thetab,psib)*ubs;
Jsnphib = Cbn(phi,theta,psi)*R1(psib)*R2(thetab)*JR3(phib)*ubs;
Jsnthetab = Cbn(phi,theta,psi)*R1(psib)*JR2(thetab)*R3(phib)*ubs;
Jsnpsib = Cbn(phi,theta,psi)*JR1(psib)*R2(thetab)*R3(phib)*ubs;
Jsngama = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*[0;cosd(gama(i));-sind(gama(i))];
Jsn = [Jsnphi Jsntheta Jsnpsi Jsnt Jsnphib Jsnthetab Jsnpsib Jsngama];
Denominateur = (S*Jrn(:,8));
a1(i) = -((S*JXnphi)/Denominateur);
a2(i) = -((S*JXntheta)/Denominateur);
a3(i) = -((S*JXnpsi)/Denominateur);
a4(i) = -((S*JXnt)/Denominateur);
a5(i) = -((S*Jrn(:,5))/Denominateur);
a6(i) = -((S*Jrn(:,6))/Denominateur);
a7(i) = -((S*Jrn(:,7))/Denominateur);
a9(i) = -((S*Jrn(:,9))/Denominateur);
coef = (a1(i)^2*sdphi^2+a2(i)^2*sdtheta^2+a3(i)^2*sdpsi^2+a4(i)^2*dTli^2+a5(i)^2*sdphib^2+a6(i)^2*sdthetab^2+a7(i)^2*sdpsib^2+a9(i)^2*sdgama^2);
vec = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*ubs;
CoVar_Sn = coef*vec*vec';
SXn_Proposed = CoVar_Pn + CoVar_an + CoVar_rn + CoVar_Sn;
SXn_Classical = CoVar_Pn + CoVar_an + CoVar_rn;
ITC_total_Proposed = sqrt(diag(SXn_Proposed));
ITC_TOTAL_Proposed(:,i) = ITC_total_Proposed;
Incertitude_3D_Proposed = sqrt(ITC_TOTAL_Proposed(1,:).^2+ITC_TOTAL_Proposed(2,:).^2+ITC_TOTAL_Proposed(3,:).^2);
ITC_total_Classical = sqrt(diag(SXn_Classical));
ITC_TOTAL_Classical(:,i) = ITC_total_Classical;
Incertitude_3D_Classical = sqrt(ITC_TOTAL_Classical(1,:).^2+ITC_TOTAL_Classical(2,:).^2+ITC_TOTAL_Classical(3,:).^2);
Comparaison_SXn = Incertitude_3D_Proposed-Incertitude_3D_Classical;
Pourcentage = Comparaison_SXn./Incertitude_3D_Proposed*100;
ITC_TOTALW(:,i) = [ITC_TOTAL_Classical(:,i);ITC_TOTAL_Proposed(:,i);Pourcentage(i)];
SXn = CoVar_Pn + CoVar_an + CoVar_rn + CoVar_Sn;
% SXn = CoVar_Pn + CoVar_an + CoVar_rn;
SXn = CoVar_Sn;

ITC_total = sqrt(diag(SXn));
ITC_TOTAL(:,i) = ITC_total;
Xn = Pn + Cbn(phi,theta,psi)*(Cbn(phib,thetab,psib)*[X;rho(i)*sind(gama(i));rho(i)*cosd(gama(i))] + abi);
Xn_TOTAL(:,i) = Xn;
CoVar_PnTOTAL(:,i) = sqrt(diag(CoVar_Pn));
CoVar_anTOTAL(:,i) = sqrt(diag(CoVar_an));
CoVar_rbsTOTAL(:,i) = sqrt(diag(CoVar_rn));
CoVar_SnTOTAL(:,i) = sqrt(diag(CoVar_Sn));
RHOi(1,i) = drho_i;
%Imprimer les ellipsoïdes d'incertitude-type composée
%  h2 = plot_gaussian_ellipsoid([Xn(1),-Xn(2),-Xn(3)],SXn,30);
 h2 = plot_gaussian_ellipsoid([-Xn(2),-Xn(3)],SXn(2:3,2:3),40);
  set(h2,'color',[1,0,1]); 


%  set(h2,'facealpha',0.6);
%  s.FaceColor = 'red';
 
 

%  grid on; axis equal; axis tight;
%  set(gcf,'color','w');
%  set(gca,'Fontsize',15);
%  hold on
%  plot3(Xn(1),-Xn(2),-Xn(3),'*');
%  plot3(Pn(1),-Pn(2),-Pn(3),'*');
%  Dist = [Pn';Xn'];
%  plot3(Dist(:,1), -Dist(:,2),-Dist(:,3),'b')
%  STOTAL(i,:) = S;

 grid on; axis equal; axis tight;
 set(gcf,'color','w');
 set(gca,'Fontsize',15);
 hold on
  plot(-Xn(2),-Xn(3),'*','Color','r');

 Dist = [Pn';Xn'];
 plot( -Dist(:,2),-Dist(:,3),'b')
 STOTAL(i,:) = S;
end
% title('Profil Longitudinal')
% XnITC = [Xn_TOTAL' ITC_TOTAL' RHOi'];
% [SV, curvature] = findPointNormals_mg(Xn_TOTAL', 10, 0.1);

Incertitude_3D = sqrt(ITC_TOTAL(1,:).^2+ITC_TOTAL(2,:).^2+ITC_TOTAL(3,:).^2);
% view(90,0); set(gca,'proj','perspective'); grid on; 
FXY = 20;
Ftitle = 22;
Posi = [-80 -60 -40 -20 0 20 40 60 80];
xticks(Posi)
xticklabels({'-80','-60','-40','-20','0','20','40','60','80'})
% set(gca,'XTickLabel',h2,'FontSize',FXY,'FontAngle', 'italic')
ax = gca;
ax.FontSize = FXY; 

Ftitle = 20;
PosiY = [-50 -40 -30 -20 -10 0];
yticks(PosiY)
yticklabels({'0','-10','-20','-30','-40','-50'})
ay = gca;
ay.FontSize = FXY; 
% xlim([-5 5])
xlim([-90 90])
ylim([-50 0])
% ylim([20 80])
grid on
% grid minor
title('Flat Seafloor','FontSize', Ftitle,'FontAngle', 'italic')
xlabel('Across-track distance (m)', 'FontSize', FXY,'FontAngle', 'italic')
ylabel('Depth (m)', 'FontSize', FXY,'FontAngle', 'italic')
hold off
% saveas(h2, 'UncertaintyFlat.png')





figure
hold on
for i=1:k
if any(isnan(Wpente(i)))
    Wpente(i)= Wpente(i-1);
else
end
alpha(i) = Wpente(i)*d2r;
%Paramètres du plan associés à l'inclinaison des dunes (alpha)
S = [0 sin(alpha(i)) cos(alpha(i))];
%**************************************************************************
xn = 0;
yn = 0;
zn = 0;
Pn = [xn;yn;zn];
%Vitesses Linéaires de la plateforme
VN = 0;
VE = 0;
VD = 0;
%Centrale inertielle (phi, theta, psi, w1, w2, w3)
%Angles d'attitude
phi = 0*d2r;
theta = 0*d2r;
psi = 0*d2r;
%Écarts-type des angles d'attitude
sdphi = 0.03*d2r;
sdtheta = 0.03*d2r;
sdpsi = 0.03*d2r;
%Vitesses Angulaires
w1 = 0*d2r;
w2 = 0*d2r;
w3 = 0*d2r;
Omega = [0 -w3 w2;w3 0 -w1;-w2 w1 0];
%--------------------------------------------------------------------------
%Matrice Variance-Covariance pour la position(Pn) -->Pn=Pn(Xn,Yn,Zn,t)
CoVar_xi = [sdxn^2 0 0 0;
         0 sdyn^2 0 0;
         0 0 sdzn^2 0;
         0 0 0 dTlp^2];
%Derivées partielles de Pn (Matrice des derivées Position)
JPn = [1 0 0 VN;
       0 1 0 VE;
       0 0 1 VD];
%Matrice Variance-Covariance Positionnement
CoVar_Pn = JPn*CoVar_xi*JPn';
%--------------------------------------------------------------------------
%Matrice Variance-Covariance lever-arms(an) -->zeta=(phi,theta,psi,t,ax,ay,az)
%Derivées partielles de an (Matrice des derivées partielles Lever-arms)
Jan = [R1(psi)*R2(theta)*JR3(phi)*abi R1(psi)*JR2(theta)*R3(phi)*abi JR1(psi)*R2(theta)*R3(phi)*abi Cbn(phi,theta,psi)*Omega*abi Cbn(phi,theta,psi)];
%MVC des paramètres de l'effet des lever-arms
CoVar_zeta = [sdphi^2 0 0 0 0 0 0;
       0 sdtheta^2 0 0 0 0 0;
       0 0 sdpsi^2 0 0 0 0;
       0 0 0 dTli^2 0 0 0;
       0 0 0 0 sdax^2 0 0;
       0 0 0 0 0 sday^2 0;
       0 0 0 0 0 0 sdaz^2];
CoVar_an = Jan*CoVar_zeta*Jan';
ubs = [0;sind(gama(i));cosd(gama(i))];
rbs = rho(i)*ubs;
Jrbs_rho =[0; sind(gama(i)); cosd(gama(i))];
Jrbs_gama =rho(i)*[0;cosd(gama(i));-sind(gama(i))];
%Écarts-types des données fournies par le MBES
wtt = (rho(i)/SVP_SIMULE)*2;
sdrho = [1/2*SVP_SIMULE 1/2*wtt]*[stdMBES^2 0;0 sCEL(1)^2]*[1/2*SVP_SIMULE; 1/2*wtt];
Jrnphi = R1(psi)*R2(theta)*JR3(phi)*Cbn(phib,thetab,psib)*rbs;
Jrntheta = R1(psi)*JR2(theta)*R3(phi)*Cbn(phib,thetab,psib)*rbs;
Jrnpsi = JR1(psi)*R2(theta)*R3(phi)*Cbn(phib,thetab,psib)*rbs;
Jrnt = Cbn(phi,theta,psi)*Omega*Cbn(phib,thetab,psib)*rbs;
Jrnphib = Cbn(phi,theta,psi)*R1(psib)*R2(thetab)*JR3(phib)*rbs;
Jrnthetab = Cbn(phi,theta,psi)*R1(psib)*JR2(thetab)*R3(phib)*rbs;
Jrnpsib = Cbn(phi,theta,psi)*JR1(psib)*R2(thetab)*R3(phib)*rbs;
Jrnrho = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*Jrbs_rho;
Jrngama = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*Jrbs_gama;
Jrn = [Jrnphi Jrntheta Jrnpsi Jrnt Jrnphib Jrnthetab Jrnpsib Jrnrho Jrngama];
%Dérivées partielles de Xn 
JXnphi = Jrn(:,1) + Jan(:,1); 
JXntheta = Jrn(:,2) + Jan(:,2); 
JXnpsi = Jrn(:,3) + Jan(:,3); 
JXnt = Jrn(:,4) + Jan(:,4); 
%Calcul de la matrice variance covariance pour rn
CoVar_chi = [sdphi^2 0 0 0 0 0 0 0 0;
       0 sdtheta^2 0 0 0 0 0 0 0;
       0 0 sdpsi^2 0 0 0 0 0 0;
       0 0 0 dTli^2 0 0 0 0 0;
       0 0 0 0 sdphib^2 0 0 0 0;
       0 0 0 0 0 sdthetab^2 0 0 0;
       0 0 0 0 0 0 sdpsib^2 0 0;
       0 0 0 0 0 0 0 sdrho^2 0;
       0 0 0 0 0 0 0 0 sdgama^2];
CoVar_rn = Jrn*CoVar_chi*Jrn';
den1 = S*Jrn(:,8);
JXn_chi = [JXnphi JXntheta JXnpsi JXnt Jrn(:,5) Jrn(:,6) Jrn(:,7) Jrn(:,8) Jrn(:,9)];
JXn_chi1 = [JXnphi JXntheta JXnpsi JXnt Jrn(:,5) Jrn(:,6) Jrn(:,7) Jrn(:,9)];
dchi1 = [sdphi; sdtheta; sdpsi; dTli; sdphib; sdthetab; sdpsib; sdgama];
num1 = S*JXn_chi1*dchi1;
drho_i =-(num1/den1);
Jsnphi = R1(psi)*R2(theta)*JR3(phi)*Cbn(phib,thetab,psib)*ubs;
Jsntheta = R1(psi)*JR2(theta)*R3(phi)*Cbn(phib,thetab,psib)*ubs;
Jsnpsi = JR1(psi)*R2(theta)*R3(phi)*Cbn(phib,thetab,psib)*ubs;
Jsnt = Cbn(phi,theta,psi)*Omega*Cbn(phib,thetab,psib)*ubs;
Jsnphib = Cbn(phi,theta,psi)*R1(psib)*R2(thetab)*JR3(phib)*ubs;
Jsnthetab = Cbn(phi,theta,psi)*R1(psib)*JR2(thetab)*R3(phib)*ubs;
Jsnpsib = Cbn(phi,theta,psi)*JR1(psib)*R2(thetab)*R3(phib)*ubs;
Jsngama = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*[0;cosd(gama(i));-sind(gama(i))];
Jsn = [Jsnphi Jsntheta Jsnpsi Jsnt Jsnphib Jsnthetab Jsnpsib Jsngama];
Denominateur = (S*Jrn(:,8));
a1(i) = -((S*JXnphi)/Denominateur);
a2(i) = -((S*JXntheta)/Denominateur);
a3(i) = -((S*JXnpsi)/Denominateur);
a4(i) = -((S*JXnt)/Denominateur);
a5(i) = -((S*Jrn(:,5))/Denominateur);
a6(i) = -((S*Jrn(:,6))/Denominateur);
a7(i) = -((S*Jrn(:,7))/Denominateur);
a9(i) = -((S*Jrn(:,9))/Denominateur);
coef = (a1(i)^2*sdphi^2+a2(i)^2*sdtheta^2+a3(i)^2*sdpsi^2+a4(i)^2*dTli^2+a5(i)^2*sdphib^2+a6(i)^2*sdthetab^2+a7(i)^2*sdpsib^2+a9(i)^2*sdgama^2);
vec = Cbn(phi,theta,psi)*Cbn(phib,thetab,psib)*ubs;
CoVar_Sn = coef*vec*vec';
SXn_Proposed = CoVar_Pn + CoVar_an + CoVar_rn + CoVar_Sn;
SXn_Classical = CoVar_Pn + CoVar_an + CoVar_rn;
ITC_total_Proposed = sqrt(diag(SXn_Proposed));
ITC_TOTAL_Proposed(:,i) = ITC_total_Proposed;
Incertitude_3D_Proposed = sqrt(ITC_TOTAL_Proposed(1,:).^2+ITC_TOTAL_Proposed(2,:).^2+ITC_TOTAL_Proposed(3,:).^2);
ITC_total_Classical = sqrt(diag(SXn_Classical));
ITC_TOTAL_Classical(:,i) = ITC_total_Classical;
Incertitude_3D_Classical = sqrt(ITC_TOTAL_Classical(1,:).^2+ITC_TOTAL_Classical(2,:).^2+ITC_TOTAL_Classical(3,:).^2);
Comparaison_SXn = Incertitude_3D_Proposed-Incertitude_3D_Classical;
Pourcentage = Comparaison_SXn./Incertitude_3D_Proposed*100;
ITC_TOTALW(:,i) = [ITC_TOTAL_Classical(:,i);ITC_TOTAL_Proposed(:,i);Pourcentage(i)];
SXn = CoVar_Pn + CoVar_an + CoVar_rn + CoVar_Sn;
% SXn = CoVar_Pn + CoVar_an + CoVar_rn;
% SXn = CoVar_Sn;

ITC_total = sqrt(diag(SXn));
ITC_TOTAL(:,i) = ITC_total;
Xn = Pn + Cbn(phi,theta,psi)*(Cbn(phib,thetab,psib)*[X;rho(i)*sind(gama(i));rho(i)*cosd(gama(i))] + abi);
Xn_TOTAL(:,i) = Xn;
CoVar_PnTOTAL(:,i) = sqrt(diag(CoVar_Pn));
CoVar_anTOTAL(:,i) = sqrt(diag(CoVar_an));
CoVar_rbsTOTAL(:,i) = sqrt(diag(CoVar_rn));
CoVar_SnTOTAL(:,i) = sqrt(diag(CoVar_Sn));
RHOi(1,i) = drho_i;
%Imprimer les ellipsoïdes d'incertitude-type composée
%  h2 = plot_gaussian_ellipsoid([Xn(1),-Xn(2),-Xn(3)],SXn,30);
%  h2 = plot_gaussian_ellipsoid([-Xn(2),-Xn(3)],SXn(2:3,2:3),30);
%   set(h2,'color',[0.6350, 0.0780, 0.1840]); 


%  set(h2,'facealpha',0.6);
%  s.FaceColor = 'red';
 
 

%  grid on; axis equal; axis tight;
%  set(gcf,'color','w');
%  set(gca,'Fontsize',15);
%  hold on
%  plot3(Xn(1),-Xn(2),-Xn(3),'*');
%  plot3(Pn(1),-Pn(2),-Pn(3),'*');
%  Dist = [Pn';Xn'];
%  plot3(Dist(:,1), -Dist(:,2),-Dist(:,3),'b')
%  STOTAL(i,:) = S;

hold on
 grid on; axis equal; axis tight;
 set(gcf,'color','w');
 set(gca,'Fontsize',15);
 hold on
  plot(-Xn(2),-Xn(3),'*','Color','r');

 Dist = [Pn';Xn'];
 plot( -Dist(:,2),-Dist(:,3),'b')
 STOTAL(i,:) = S;
end
FXY = 20;
Ftitle = 22;
Posi = [-80 -60 -40 -20 0 20 40 60 80];
xticks(Posi)
xticklabels({'-80','-60','-40','-20','0','20','40','60','80'})
% set(gca,'XTickLabel',h2,'FontSize',FXY,'FontAngle', 'italic')
ax = gca;
ax.FontSize = FXY; 

Ftitle = 20;
PosiY = [-50 -40 -30 -20 -10 0];
yticks(PosiY)
yticklabels({'0','-10','-20','-30','-40','-50'})
ay = gca;
ay.FontSize = FXY; 
% xlim([-5 5])
xlim([-90 90])
ylim([-50 0])
% ylim([20 80])
grid on
% grid minor
title('Flat Seafloor','FontSize', Ftitle,'FontAngle', 'italic')
xlabel('Across-track distance (m)', 'FontSize', FXY,'FontAngle', 'italic')
ylabel('Depth (m)', 'FontSize', FXY,'FontAngle', 'italic')
hold off