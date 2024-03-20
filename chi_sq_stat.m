function [v, std]=chi_sq_stat(er, vel, dof)

% Function calculates velociti estimates and its standard deviation using
% chi squared statistics
%
% INPUTS:
%    er      - chi squared field as vector with N elements
%    vel     - velocity field used to calculate chisquared as (Nx2) matrix
%    dof     - degrees of freedom
%
% OUTPUTS:
%    v    - 2D velocity estimate
%    std  - 2D velocity standard deviation
%

  er=er/min(er);
  er=er*(dof);
  er=er-(dof);

  vv = vel(min(er)==er,:);
  v = vv(1,:);

  V=vel;
  E=er;

  vel=V(14:end,:);
  er=E(14:end);

  warning('off','all')

  [xGrid, yGrid] = meshgrid(min(vel(:,1)):0.01:max(vel(:,1)),...
      min(vel(:,2)):0.01:max(vel(:,2)));

  F = scatteredInterpolant(vel(:,1),vel(:,2),er');

  zGrid = F(xGrid,yGrid);

  I=find(zGrid<=2.3);
  
  if length(I)==1
      zmin = min(er(er~=min(er)));
      I=find(zGrid<=zmin);
  end

  Sx = max(abs([v(1)-min(xGrid(I)) max(xGrid(I))-v(1)]));
  Sy = max(abs([v(2)-min(yGrid(I)) max(yGrid(I))-v(2)]));

  ssx = max([abs(min(vel(:,1))-v(1)) abs(v(1)-max(vel(:,1)))]);
  ssy = max([abs(min(vel(:,2))-v(2)) abs(v(2)-max(vel(:,2)))]);

  if isempty(Sx) || isempty(Sy)
      Sx = ssx;
      Sy = ssy;
  end

  if Sx == ssx || Sy == ssy

      [xGrid, yGrid] = meshgrid(v(1)-2:0.01:v(1)+2,...
          v(2)-2:0.01:v(2)+2);
      F = scatteredInterpolant(V(:,1),V(:,2),E');
      zGrid = F(xGrid,yGrid);

      I=find(zGrid<=1);

      Sx = max(abs([v(1)-min(xGrid(I)) max(xGrid(I))-v(1)]));
      Sy = max(abs([v(2)-min(yGrid(I)) max(yGrid(I))-v(2)]));
  end

  if isempty(Sx) || isempty(Sy)
      Sx = ssx;
      Sy = ssy;
  end

  std = [Sx Sy];

  warning('on','all')

end