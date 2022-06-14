function u = divdiff2(x, y, n)
%!!
%!! Computes divided difference using the expanded form
%!!

  n = length(x);

  u = 0.0;
  for j=1:n
     tmp = 1.0;
     for k=1:n
       if(k ~= j) 
         tmp = tmp * (x(j)-x(k));
       end
     end
     u = u + y(j)/tmp;
  end

end % function divdiff2



