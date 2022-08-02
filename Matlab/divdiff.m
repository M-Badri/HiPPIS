
function table =  divdiff(x, y, d)
%! This subroutine computes the table of divided differences
% 
% INPUT:
% x: 1D vector of x coorrdinates
% y: 1D vector of y coorrdinates. y are the data values associated to the locations x.
% d: maximum number of consecutive mesh points used to compute divided differences.
% n: number of points in x.
%
% OUTPUT:
%
% table: array of dimension of n X (d+1)
% table = u[1], u[1,2], u[1,3], ... u[1,d+1]
%         u[2], u[2,3], u[2,4], ... u[2,d+2]
%         u[3], u[3,4], u[3,5], ... u[2,d+3]
%          .  ,    .  ,    .  , ...    .
%          .  ,    .  ,    .  , ...    .
%          .  ,    .  ,    .  , ...    .
%         u[n],    0  ,    0  , ...    .

  n = length(x);
  table = zeros(n, d+1);


  for i=1:n
    table(i,1) = y(i);
  end

  for j=2:d+1
    for i=1:n-(j-1)
      %tmp = divdiff2(x(i:i+j-1), y(i:i+j-1), j);
      %table(i,j) = tmp;

      table(i,j) = (table(i+1, j-1)-table(i, j-1)) / (x(i+j-1)-x(i));
      %!!if(abs(table(i, j) - tmp) > 1e-10 ) then
      %!!  write(*,*) 'TAJO tmp =', tmp, 'table(i,j) =', table(i,j)
      %!!  write(*,*) 'i=',i, 'j=', j
      %!!endif
    end
  end

end % function 
 

