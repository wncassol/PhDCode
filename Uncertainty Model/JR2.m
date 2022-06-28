function [JR2_theta] =JR2(theta)
JR2_theta = [-sin(theta) 0 cos(theta);0 0 0;-cos(theta) 0 -sin(theta)];