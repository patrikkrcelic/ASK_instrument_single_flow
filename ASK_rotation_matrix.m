function [rot, x0, y0, z0] = ASK_rotation_matrix(El, Az, zenith, siz)

% Function creates rotation matrix needed to transfer coordinates between
% geographic coordinates and tomography reference frame. Additionaly, 
% function calculates the geographic coordinates of the center of the 
% tomography volume. 
%
% INPUTS:
%    El     - Matrix of pixel elevation values
%    Az     - Matrix of pixel azimuth values
%    zenith - pixel position of zenith, 2D vector
%    siz     - size of the image
%
% OUTPUTS:
%    rot - rotation matrix
%    x0  - x coordinate of the center of the tomography volume
%    y0  - y coordinate of the center of the tomography volume
%    z0  - z coordinate of the center of the tomography volume
%


  th=pi/2-El(zenith(1), zenith(2));
  lam=Az(zenith(1), zenith(2));
    
  rot1=[cos(-lam) -sin(-lam) 0; sin(-lam) cos(-lam) 0; 0 0 1];
  rot=[1 0 0; 0 cos(-th) -sin(-th); 0 sin(-th) cos(-th)];

  rot=rot1*rot;

  z0=210000;
  x0=z0*tan(pi/2-El(round(siz/2),round(siz/2)))*sin(Az(round(siz/2),round(siz/2)));
  y0=z0*tan(pi/2-El(round(siz/2),round(siz/2)))*cos(Az(round(siz/2),round(siz/2)));

  Rt=[x0 y0 z0]*rot;
  x0=Rt(1);
  y0=Rt(2);
  z0=Rt(3);

end