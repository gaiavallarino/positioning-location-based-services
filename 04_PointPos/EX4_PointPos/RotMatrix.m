function [RM] = RotMatrix(phi,lambda)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
RM = [-sin(lambda) cos(lambda) 0; -sin(phi)*cos(lambda) -sin(phi)*sin(lambda) cos(phi); cos(phi)*cos(lambda) cos(phi)*sin(lambda) sin(phi)]; 

end

