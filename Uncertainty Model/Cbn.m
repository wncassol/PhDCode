function [Cbn] =Cbn(phi,theta,psi)
Cbn = R1(psi)*R2(theta)*R3(phi);