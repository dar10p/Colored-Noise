function [ x ] = f_alpha_uniform( n, range, alpha )

% Purpose:
%     Generates a discrete colored noise vector of size n with power 
%     spectrum distribution of alpha
%     White noise is sampled from Uniform (-range,range) distribution
%
% Usage:
%        [ x ] = f_alpha_uniform( n, range, alpha )
%     
%     n - problem size
%     range - range of the underlying Uniform distribution
%     alpha - resulting colored noise has 1/f^alpha power spectrum



%
%  Generate the coefficients Hk.
%

  hfa = zeros ( 2 * n, 1 );
  hfa(1) = 1.0; 
  for i = 2 : n
    hfa(i) = hfa(i-1) * ( 0.5 * alpha + ( i - 2 ) ) / ( i - 1 );
  end
  hfa(n+1:2*n) = 0.0;
  
%
%  Fill Wk with white noise.
%
  
  wfa = [ -range + 2 * range.* rand( n, 1 ); zeros( n, 1 ); ];
  
%
%  Perform the discrete Fourier transforms of Hk and Wk.
%

  [ fh ] = fft( hfa );
  [ fw ] = fft( wfa );
  
%
%  Multiply the two complex vectors.
%
    fh = fh( 1:n + 1 );
    fw = fw( 1:n + 1 );
      
    fw = fh .* fw;
    
%
%  This scaling is introduced only to match the behavior
%  of the Numerical Recipes code...
%

    fw(1)       = fw(1) / 2;
    fw(end)     = fw(end) / 2;
    
%
%  Take the inverse Fourier transform of the result.
%
  
  fw = [ fw; zeros(n-1,1); ];
  
  x = ifft( fw );
  
  x = 2*real( x(1:n) );
  
%
%  Discard the second half of the inverse Fourier transform.
%

  return
end