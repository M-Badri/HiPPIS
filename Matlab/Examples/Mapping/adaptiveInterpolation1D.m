function [yout, deg, uumin, uumax, xtable, u_table, lambda_table, sigma_table, prod_sigma_table] = adaptiveinterpolation1D(x, y, xout, degree, limiter, ... 
                                         uumin, uumax, umin_umax_input, eps0)
%
%This function adaptively build interpolants that are the evaluated at xout.
%
% INPUT: 
% n: the number points in the 1D vector x.
% m: the number of points in the 1D vector xout.
% x: 1D mesh points of length n. For i=1, ..., n-1 x_{i} >  x_{i+1}
% y: 1D vector that have the data values associated with the points x_{i} for i=1, ..., n
% xout: 1D vector of length m that represent the locations where we which the interpolate to.
% limiter: used to determine the type limiter to be used to build interpolant.
%   - limiter=0: the interpolants are built by adaptively selecting the stencil points. This a essentially 
%                non-oscillatory (ENO) approach.[https://www.sciencedirect.com/science/article/pii/S0021999196956326]).
%   - limiter=1: a data-bounded limiter is used to ensure that the interpolant is bounded by the data values.
%                The data bounded limiter is base on Berzins' work [https://epubs.siam.org/doi/pdf/10.1137/050625667].
%   - limiter=2: a positivity-preserving limiter is used ensure that interpolant is not bounded by the data values but
%                remains positive.
% degree: degree of desired interpolant for each interval.
%
% OUTPUT:
% yout: results of evaluating interpolants at the locations xout.
% deg: 1D vector that holds the degree of the interpolant used for each interval
%

  n = length(x);
  m = length(xout);
  yout = zeros(size(xout));
  u_table = zeros(n-1, degree+1);           % table of devided diferences
  lambda_table = zeros(n-1, degree+1);      % table of lambda vlaues
  sigma_table = zeros(n-1, degree+1);       % table of sigma values    
  prod_sigma_table = zeros(n-1, degree+1);    %       
  x_table = zeros(n-1, degree+1);           % table of stencil choice 
  slope = zeros(n-1, 1);
  mm_l = zeros(n-1, 1);
  mm_r = zeros(n-1, 1);
  www = zeros(n-1, 1);
  u = zeros(degree+1, 1);
  xval =zeros(degree+1, 1); 
  up_b = zeros(degree+1, 1);
  low_b = zeros(degree+1, 1);
  prod_deltax = zeros(degree+1, 1);
  lambda = zeros(degree+1, 1);
  sigma = zeros(degree+1, 1);
  prod_sigma = zeros(degree+1, 1);
  e = zeros(degree+1, 1);
  m_sigma = 100.0;




  %!!** Check input to make sure x_{i} < x_{i+1} **!!
  for i=1:n-1
    if( x(i) >= x(i+1) || abs(x(i+1)-x(i)) <= 1e-16 )
     fprintf('ERROR: Incorrect input at i=, %d, x(i)= %d, x(i+1)= %d \n',i, x(i), x(i+1));
     fprintf('x(i) must be less that x(i+1) and |x(i+1)-x(i)| mus be greater than machine precision eps.' );
     stop
     %exit
    end
  end

  %!!** Check output to make sure x_{i} < x_{i+1} **!!
  for i=1:m-1
    if( xout(i) >= xout(i+1) ) 
      fprintf('ERROR: Incorrect output at k= %d, xout(k) = %d, xout(k+1)= %d must be less that xout(k+1) \n', i, xout(i), xout(i+1) );
      stop
      %exit
    end
  end


  %!!** Initialize variables **!!
  umin_umax_type = 2;            %% used to select how to compute bounds u_{min} and u_{max} 
  stencil_type = 3;              %% used to choose the stencil selection process
  tol = 1.0e-8;
  table = 0.0;  
  eps = 1e-16;                   %% defined as epsilon 
  inv_eps = 1e+16;               %% defined to be + infinity
  k = 1;                         %% iteration idex used for output points

  %!!** Using eps0 to set eps2. eps0 is a user defined parameter used in the PPI
  %!!   method to relax the bounds the interpolant for the cases where hiden
  %!!   extremum is dectected. eps0 and eps2 should be set to small values
  %!!   otherwise this may lead to large oscillation for the PPI algorithm **!!
  if (exist('eps0'))    
    eps2 = eps0;
  else
    eps2 = 0.010;
  end
  if(exist('x_table') )
    x_table = 0.0;
  end
  if(exist('u_table') )
    u_table = 0.0;
  end
  if(exist('lambda_table') )
    lambda_table = 0.0;
  end
  

  table = divdiff(x, y, degree+3);    %compute the table of divided differences 


  %!!** Calculate slopes for each interval **!!
  slope(1) = (y(3)-y(2))/(x(3)-x(2));  %!! left boundary
  for i=1:n-1
    slope(i+1) = (y(i+1)-y(i))/(x(i+1)-x(i));  %!! right boundary
  end
  slope(n+1) = slope(n-1);
  

  %!!** Mark intervals with extrema **!!
  %extrema(1) = 0
  %for i=1, n-1
  %  if(slope(i+1)*slope(i) < 0 .or. slope(i)*slope(i+2) < 0 ) then
  %    extrema(i+1) = 1   
  %  else
  %    extrema(i) = 0
  %  end
  %end
  %extrema(n+1) = 0

  %!!** Calculate the polynomial bounds for each interval **!!
  if(degree > 1 && limiter == 2)
  for i=1:n-1
    %!!** the slople for interval i **!!
    slope_im1 = slope(i);
    slope_i   = slope(i+1);
    slope_ip1 = slope(i+2);

    umin = min(y(i), y(i+1));
    umax = max(y(i), y(i+1));
    u2= (y(i+1)-y(i))/(x(i+1)-x(i));    %!! set second slected divided difference

    %!!** Option 1: umin_umax_type == 1. Note that if y(i) and y(i+1) are positive 
    %!!   u_min = 0.0 and u_max = 2* max(y(i), y(i+1)). **!!
    if(umin_umax_type == 1)
      tmp1 = min(y(i), y(i+1));
      tmp2 = max(y(i), y(i+1));
      umin = tmp1 - abs(tmp1);
      umax = tmp2 + abs(tmp2);
    end
  
    %!!** Option 2: umin_umax_type == 2. In intervals where a hidden local 
    %!!   extremum is detectected  umin and/or umax are Calculation of umin umax option 2 **!!
    if(umin_umax_type == 2) 

      tmp1 = min(y(i), y(i+1));
      tmp2 = max(y(i), y(i+1));

      %!!** Calculcualtion of umin based of the existence of an extremum **!!
      if( (slope_im1*slope_ip1 < 0.0 && slope_im1 < 0.0) ||     ...      %% Detects a minimum
          (slope_im1*slope_ip1 > 0.0 && slope_im1*slope_i < 0.0) )  %% Detects a maximum and/or minimum (ambiguous).
        umin = tmp1 - abs(tmp1);
      else                                                                  % No extremum detected
        umin = tmp1 - eps2*abs(tmp1);
      end 

      %!!** Calculation of umax based on the existence of an extremum **!!
      if( (slope_im1*slope_ip1 < 0.0 && slope_im1 > 0.0)  ||  ...          % Detects a maximum
          (slope_im1*slope_ip1 > 0.0 && slope_im1*slope_i < 0.0) )   % Detects a minimum and/or a maxmimum
        umax = tmp2 + abs(tmp2);
      else                                                                  % No extremum is detected
        umax = tmp2 + eps2*abs(tmp2);
      end

    end
    
    if (exist('uumin') )
      if( exist('umin_umax_input') ) 
        %!!** In case uumin and uumax are provided, used use provided values to
        %!!   update umin and umax **!! 
        if( umin_umax_input)
          umin = uumin(i);
        else
        %!!** Incase uumin and uumax are not provide but the arrays uumin and uumax
        %!!   are passed update the arrays using umin and umax **!! 
          uumin(i) = umin;
        end
      end
    end
    if (exist('uumax') )
      if( exist('umin_umax_input') )
        %!!** In case uumin and uumax are provided, used use provided values to
        %!!   update umin and umax **!! 
        if( umin_umax_input ) 
          umax = uumax(i);
        else
        %!!** Incase uumin and uumax are not provide but the arrays uumin and uumax
        %!!   are passed update the arrays using umin and umax **!! 
          uumax(i) = umax;
        end
      end
    end

    %%!!** Compute the values of m_{\ell} and m_r for the positive-preserving 
    %%!!   method. This coresponds to the default setting **!!
    if(y(i) < y(i+1)) 
      ww = u2;
      m_l = (umin-y(i)) / (y(i+1)-y(i));
      m_l = min(0.0, m_l);
      m_r = (umax-y(i)) /  (y(i+1)-y(i)); 
      m_r = max(1.0, m_r);
    elseif(y(i) > y(i+1))
      ww = u2;
      m_l = (umax-y(i)) /  (y(i+1)-y(i)); 
      m_l = min(0.0, m_l);
      m_r = (umin-y(i)) / (y(i+1)-y(i));  
      m_r = max(1.0, m_r);
    else % This part deals with the special case where y(i)=y(i+1) 
      ww = u2;
      tmp_si = max(i-1, 1);
      ul = table(tmp_si, 3 );    % divided difference U[x_{i-1}, x_{i}, x_{i+1}] 
      tmp_ei = min(i+2,n); 
      ur = table(i, min(i, tmp_ei-i+1) );  % divided difference U[x_{i}, x_{i+1}, x_{i+2}] 
      if( ul >  0.0 && ur ~= 0.0)
        ww = ul * (x(i+1)-x(i)) * (x(i+1)-x(tmp_si)); 
        m_l = (umin-y(i)) / ww ; 
        m_l = min(0.0, m_l);
        m_r = (umax-y(i)) / ww; 
        m_r = max(1.0, m_r);
      elseif (ul < 0.0  && ur ~= 0.0)
        ww = ul * (x(i+1)-x(i)) * (x(i+1)-x(tmp_si)) ;
        %u2 = ww
        m_l = (umax-y(i)) /  ww ; 
        m_l = min(0.0, m_l);
        m_r = (umin-y(i)) / ww;  
        m_r = max(1.0, m_r);
      elseif (ur > 0.0 && ul ~= 0.0) 
        ww = ur * (x(i+1)-x(i)) * (x(tmp_ei)-x(i)) 
        %u2 = ww
        m_l = (umin-y(i)) / ww; 
        m_l = min(0.0, m_l);
        m_r = (umax-y(i)) / ww; 
        m_r = max(1.0, m_r);
      elseif (ur < 0.0 && ul ~= 0.0)
        ww = ur * (x(i+1)-x(i)) * (x(tmp_ei)-x(i)); 
        !!u2 = ww;
        m_l = (umax-y(i)) /  ww ; 
        m_l = min(0.0, m_l);
        m_r = (umin-y(i)) / ww;  
        m_r = max(1.0, m_r);
      else % ur=0 or ul=0, the algorithm defaults to BDI.  
        ww = u2;
        m_l = 0.0;
        m_r = 1.0;
      end
    end

    %!!** Save the bunds for each interval and the www ***!!
    mm_l(i) = m_l;
    mm_r(i) = m_r;
    www(i) = ww;
  end
  %!!** Default case: DBI **!!
  else  
    for i=1:n-1
      %!!** Compute the values of m_{\ell} and m_r for the data-bounded method 
      %!!   if the limiter variable is set to 1 **!!
      mm_l(i) = 0.0;
      mm_r(i) = 1.0;
      www(i) = (y(i+1)-y(i)) / (x(i+1)-x(i));
    end
  end

  %!!** loop over each input intervals. For each  interval build an interpolant and 
  %!!   evaluate the interpolant at the desired output points **!!
  for i=1:n-1
 
    %!!** Initialize varibles for each interval**!!
    u = zeros(degree+1,1);
    xval = zeros(degree+1,1);
    up_b = zeros(degree+1,1);
    low_b = zeros(degree+1,1);
    prod_deltax = zeros(degree+1,1);
    lambda = zeros(degree+1,1);
    sigma = ones(degree+1,1);
    prod_sigma = zeros(degree+1,1);
   
    xval(1) = x(i);                       % first point in stencil
    xval(2) = x(i+1);                     % second point in stencil
    u(1) = y(i);                          % set first selected divided difference
    u(2)= (y(i+1)-y(i))/(x(i+1)-x(i));    % set second slected divided difference
    lambda(1) = 1.0;                      % set first ratio of divided difference  
    lambda(2) = 1.0;                      % set second ratio of divided difference  
    sigma(1) = 1.0;
    sigma(2) = u(2)*(x(i+1)-x(i));
    prod_sigma(1) = 1.0;
    prod_sigma(2) = (x(i+1)-x(i));
    up_b(1) = 1.0;
    up_b(2) = 1.0;
    low_b(1) = -1.0;
    low_b(2) = -1.0;
    m_l = mm_l(i);
    m_r = mm_r(i);
    ww = www(i);
    ei=i+1;                               % set e before entering the do loop
    si=i;                                 % set s before entering the do loop
   
    %!!** Continue to build high degree interpolant only if the target degree is
    %!!   greater than 1 and the both U{x_{i}, x_{i+1}, x_{i+2}} and U{x_{i-1},x_{i}, x_{i+1}}
    %!!   are not zeros. **!!
    if(degree > 1 && ww ~= 0.0)
      for j=2:degree 

         %!!** Initialize selection boolean variables **!!
         bool_left = false;
         bool_right = false;

         tmp_si = max(1, si-1);                          % decrementing stencil left idex
         tmp_ei = min(n, ei+1);                          % incrementing stencil right index
         
         %!!** Calculate ul and ur (left and right divided differences
         %!!   respectively)
         %!!   Calculate lambda_l and lambda_r ( left and right ratio of 
         %!!   of divided references respectively ) **!!
         prod_deltax_l = prod_deltax(j) * (x(ei)-x(tmp_si)); % product of interval containing stencil
         if(xval(j) < x(i))
           prod_sigma_l = prod_sigma(j)*(x(i+1)-xval(j));
         else
           prod_sigma_l = prod_sigma(j)*(x(i)-xval(j));
         end
         prod_sigma_l2 = prod_sigma_l*(x(i+1)-x(tmp_si));
         prod_deltax_l = prod_deltax(j) * (x(ei)-x(tmp_si)); % product of interval containing stencil
         if(si-1 > 0)
           ul = table(tmp_si, ei-tmp_si+1);   % get left divided difference   
           lambda_l = ul/ww * prod_deltax_l;  % calculate left lambda         
           sigma_l = ul * prod_sigma_l;
           %sigma_l =  prod_sigma_l;
         else
           ul = inv_eps;  % set left dividided difference to + infinity  
           lambda_l = inv_eps;    % set left lambda to infinity
           sigma_l = inv_eps;
         end

         prod_deltax_r = prod_deltax(j) * (x(tmp_ei)-x(si));    % product of inteval containing stencil
         if(xval(j) < x(i) )
           prod_sigma_r = prod_sigma(j)*(x(i+1)-xval(j));
         else
           prod_sigma_r = prod_sigma(j)*(x(i)-xval(j));
         end
     
         prod_sigma_r2 = prod_sigma_r*(x(tmp_ei)-x(i));
         if(ei+1 <= n)
           ur = table(si, tmp_ei-si+1);         % get right divided difference 
           lambda_r =  ur/ ww * prod_deltax_r;  % calculate righ lambda
           sigma_r = ur * prod_sigma_r;
           %sigma_r = prod_sigma_r
         else
           ur = inv_eps;          % set righl divided difference to infinity
           lambda_r =  inv_eps;   % set righ lambda to infinity
           sigma_r = inv_eps;
         end

         e(j) = -(x(i)-xval(j))/(x(i+1)-x(i));   %calculate e_j !!
         d_l = (x(ei)-x(tmp_si))/(x(i+1)-x(i));  %calculate d_l 
         d_r = (x(tmp_ei)-x(si))/(x(i+1)-x(i));  %calculate d_r

         %!!** In the case where the points inserted to 
         %!!   V_{j-1} form V_{j} is to the left. Calculate 
         %!!   upper and lower bounds up_b_l and low_bl 
         %!!   respectively  **!!
         if(j == 2 ) 
           up_b_l = d_l*( -m_l*4.0 + 1.0 );
           low_b_l = d_l*( -(m_r-1.0)*4.0 - 1.0 ) ;
         else
           if(e(j) <= 0.0)
             up_b_l = (up_b(j) - lambda(j))* d_l / (1.0-e(j));
             low_b_l = (low_b(j) - lambda(j))*d_l / (1.0-e(j));
           elseif(e(j) > 0.0)
             up_b_l = (low_b(j) - lambda(j))*d_l / (0.0-e(j));
             low_b_l = (up_b(j) - lambda(j))* d_l / (0.0-e(j));
           end
         end
     
         %!!** In the case where the points inserted to 
         %!!   V_{j-1} form V_{j} is to the right. Calculate 
         %!!   upper and lower bounds up_b_r and low_b_r 
         %!!   respectively  **!!
         if(j == 2)
            up_b_r = d_r*( -m_l*4.0 + 1.0 );
            low_b_r = d_r*( -(m_r-1.0)*4.0 - 1.0 );  
         else
           if(e(j) <= 0.0) 
             up_b_r = (up_b(j) - lambda(j))* d_r / (1.0-e(j));
             low_b_r = (low_b(j) - lambda(j))*d_r / (1.0-e(j));
           elseif(e(j) > 0.0) 
             up_b_r = (low_b(j) - lambda(j))*d_r / (0.0-e(j));
             low_b_r = (up_b(j) - lambda(j))* d_r / (0.0-e(j));
           end
         end
        

         %!!** Option 1: stencil_type = 1. In addition to positivity or 
         %!!   data boundedness, the stencil selection is based on the ENO approach **!!
         if(stencil_type == 1) 
           if( (low_b_l <= lambda_l && lambda_l <= up_b_l) && ... %% Adding a point to left meets the requiremenst for DBI or PPI
               (low_b_r <= lambda_r && lambda_r <= up_b_r) )      %% Adding a point to right meets the requiremenst for DBI or PPI
             %!!** boolean variable is set to true based on the coresponding 
             %!!   divided difference |ul| or |ur| is the smalest
             if(abs(ul) < abs(ur) )
               bool_left = true;
               bool_right = false;
             else
               bool_left = false;
               bool_right = true;
             end
           elseif(low_b_r <= lambda_r && lambda_r <= up_b_r) % Adding a point to right meets the requiremenst for DBI or PPI
             bool_left = false;
             bool_right = true;
           elseif(low_b_l <= lambda_l && lambda_l <= up_b_l) % Adding a point to left meets the requiremenst for DBI or PPI
             bool_left = true;
             bool_right = false;
           end
         end 

         %!!** Option 2: stencil_type = 2. In adddition to DBI or PPI the 
         %!!   stencil selection is based on the smalest lambda **!!
         if(stencil_type == 2) 
           if( (low_b_l <= lambda_l && lambda_l <= up_b_l)  &&  ...
               (low_b_r <= lambda_r && lambda_r <= up_b_r) )
             if(abs(lambda_l) < abs(lambda_r) )
               bool_left = true;
               bool_right = false;
             else
               bool_left = false;
               bool_right = true;
             end
           elseif(low_b_r <= lambda_r && lambda_r <= up_b_r)  
             bool_left = false;
             bool_right = true;
           elseif(low_b_l <= lambda_l && lambda_l <= up_b_l)  
             bool_left = true;
             bool_right = false;
           end
         end 

     
         %!! Option 3: stencil_type = 3. In addition to DBI or PPI the 
         %!! stencil selection prioritize a symetric stencil other others **!!
         if(stencil_type == 3) 
           if( (low_b_l <= lambda_l && lambda_l <= up_b_l)  &&  ...
               (low_b_r <= lambda_r && lambda_r <= up_b_r) )
             if( i-si < ei-i )
               bool_left = true;
               bool_right = false;
             elseif( i-si > ei-i ) 
               bool_left = false;
               bool_right = true;
             else
               if(abs(lambda_l) < abs(lambda_r) )
                 bool_left = true;
                 bool_right = false;
               else
                 bool_left = false;
                 bool_right = true;
               end
             end
           elseif(low_b_r <= lambda_r && lambda_r <= up_b_r) 
             bool_left = false;
             bool_right = true;
           elseif(low_b_l <= lambda_l && lambda_l <= up_b_l) 
             bool_left = true;
             bool_right = false;
           end
         end 

         %!!** Stencil choice option 4 **!!
         if(stencil_type == 4) 
           if( (low_b_l <= lambda_l && lambda_l <= up_b_l)  &&  ...
               (low_b_r <= lambda_r && lambda_r <= up_b_r) )
             if( abs(sigma_l) < abs(sigma_r) )
               bool_left = true;
               bool_right = false;
             elseif( abs(sigma_l) >= abs(sigma_r) ) 
               bool_left = false;
               bool_right = true;
             else
               if(abs(lambda_l) < abs(lambda_r) )
                 bool_left = true;
                 bool_right = false;
               else
                 bool_left = false;
                 bool_right = true;
               end
             end
           elseif(low_b_r <= lambda_r && lambda_r <= up_b_r) 
             bool_left = false;
             bool_right = true;
           elseif(low_b_l <= lambda_l && lambda_l <= up_b_l) 
             bool_left = true;
             bool_right = false;
           end
         end


         if(stencil_type == 5) 
           if( (low_b_l <= lambda_l && lambda_l <= up_b_l)  &&  ...
               (low_b_r <= lambda_r && lambda_r <= up_b_r) )
             if( abs(prod_sigma_l2) < abs(prod_sigma_r2) )
               if(abs(sigma_l) < m_sigma*abs(sigma(j)) )
                 bool_left = true;
                 bool_right = false;
               elseif( abs(sigma_r) < m_sigma*abs(sigma(j)) )
                 bool_left = false;
                 bool_right = true;
               else
                 bool_left = true;
                 bool_right = false;
               end
             elseif( abs(prod_sigma_l2) > abs(prod_sigma_r2) ) 
               if( abs(sigma_r) < m_sigma*abs(sigma(j)) )
                 bool_left = false;
                 bool_right = true;
               elseif(abs(sigma_l) < m_sigma*abs(sigma(j)) ) 
                 bool_left = true;
                 bool_right = false;
               else
                 bool_left = false;
                 bool_right = true;
               end
             else
               if(abs(sigma_l) < abs(sigma_r) )
                 bool_left = true;
                 bool_right = false;
               else
                 bool_left = false;
                 bool_right = true;
               end
             end
           elseif(low_b_r <= lambda_r && lambda_r <= up_b_r) 
             bool_left = false;
             bool_right = true;
           elseif(low_b_l <= lambda_l && lambda_l <= up_b_l) 
             bool_left = true;
             bool_right = false;
           end
         end


         %!!** Add point to the left of current stencil and corresponding 
         %!!   variables **!!
         if( (bool_left == true) && (bool_right == false)) 
           si = max(1, si-1);
           %ei = ei;
           lambda(j+1) = lambda_l;
           u(j+1) = ul;
           xval(j+1) = x(si);
           up_b(j+1) = up_b_l;
           low_b(j+1) = low_b_l;
           prod_deltax(j+1) = prod_deltax_l;
           sigma(j+1) = sigma_l;
           prod_sigma(j+1) = prod_sigma_l;
         %!!** Add point to the right of current stencil and corresponding 
         %!!   variables **!!
         elseif( (bool_left == false) && (bool_right && true) ) 
           %si = si
           ei = min(ei+1, n);
           lambda(j+1) = lambda_r;
           u(j+1) = ur;
           xval(j+1) = x(ei);
           up_b(j+1) = up_b_r;
           low_b(j+1) = low_b_r;
           prod_deltax(j+1) = prod_deltax_r;
           sigma(j+1) = sigma_r;
           prod_sigma(j+1) = prod_sigma_r;
         end
         
      end % j loop
    end % end of if j > 1


    %--%!!** save the interpolant degree used for the interval [x_{i}, x_{i+1}] **!!
    deg(i) = ei-si;


    %!!!!!!** **!!
    %!!if(degree > 4)then
    %!!do j=4, degree
    %!!  if(abs(sigma(j-2)) < abs(sigma(j-1)) .and. abs(sigma(j-1)) < abs(sigma(j)) .and. abs(sigma(j)) < abs(sigma(j+1)) )then
    %!!  !!if(abs(sigma(j-1)) < abs(sigma(j)) .and. abs(sigma(j)) < abs(sigma(j+1)) )then
    %!!  !!if(abs(sigma(j)) > 10.0*abs(sigma(j+1)) ) then
    %!!    u(j-2:degree+1) = 0.0
    %!!    deg(i) = j-2
    %!!    exit
    %!!  end
    %!!end
    %!!end

    %!!** Extrapolate to points that are to the left of the defined interval **!! 
    if(k <=m)
      while( xout(k) < x(1) )
        fprintf('WARNING: Some of the output are obtained via extrapolation instead of interpolation \n.')
        fprintf(' The desired proprety such as data-boundedness or positvity is not preserved in such case \n')
        fprintf( 'k=%d, 1, x(1)= %d, xout(k) =%d \n', k, x(1), xout(k) ) ;
        yout(k) = newtonPolyVal(xval, u, xout(k));
        k = k+1;
        if(k > m) break; end
      end
    end
 
    %!!** Building and evaluating Interpolant at xout points **!!
    %!! - do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) .and. k <= m )
    if( k <=m)
      while( x(i) <= xout(k) && xout(k) <= x(i+1) )
        yout(k) = newtonPolyVal(xval, u, xout(k));
        k = k+1;
        if(k > m) break; end
      end
    end

    %!!** Extrapolate to points that are to the right of the defined interval **!! 
    if(k <= m)
      while( xout(k) > x(n) )
        fprintf('WARNING: Some of the output are obtained via extrapolation instead of interpolation. \n')
        fprintf('The desired proprety such as data-boundedness or positvity is not preserved in such case \n')
        fprintf( 'k=%d, n, x(n)= %.8E, xout(k) =%.8E \n', k, x(n), x(n)-xout(k) ); 
          yout(k) = newtonPolyVal(xval, u, xout(k));
          k = k+1;
          if(k > m) break; end
      end
    end

    %!!** Save data into table for study **!!
    if (exist('x_table') )
    for kk=1:degree+1
      x_table(i, kk) = xval(kk);
    end
    end
    if (exist('u_table') )     
    for kk=1:degree+1
      u_table(i, kk) = u(kk);
    end
    end
    if (exist('lambda_table') )      
    for kk=1:degree+1
      lambda_table(i, kk) = lambda(kk);
    end
    end
    if (exist('sigma_table') )     
    for kk=1:degree+1
      sigma_table(i, kk) = sigma(kk);
    end
    end
    if (exist('prod_sigma_table') )     
    for kk=1:degree+1
      prod_sigma_table(i, kk) = prod_sigma(kk);
    end
    end

  end %!! of i loop

end %adaptiveInterpolation1D

