function [ x ] = f_bialpha_gaussian( n, q_d, alpha1,alpha2,l0)

% alpha2 and l0 are not used in the code
% confirmed by Leandro
% Purpose:
%     Generates a discrete colored noise vector of size n with power 
%     spectrum distribution of alpha1
%     White noise is sampled from Gaussian (0,q_d) distribution
%
% Usage:
%        [ x ] = f_alpha1_gaussian( n, q_d, alpha1 )
%     
%     n - problem size
%     q_d - variance of the underlying zero-mean Gaussian
%     alpha1 - resulting colored noise has 1/f^alpha1 power spectrum


%  Set the deviation of the noise.
%
  q_d = sqrt ( q_d );
%
%  Generate the coefficients Hk.
%

  hfa = zeros ( 2 * n, 1 );
  hfa(1) = 1.0; 
  for i = 2 : n
    hfa(i) = hfa(i-1) * ( 0.5 * alpha1 + ( i - 2 ) ) / ( i - 1 );
  end
  hfa(n+1:2*n) = 0.0;
  
%
%  Fill Wk with white noise.
%
  
   wfa = [ q_d * randn( n, 1 ); zeros( n, 1 ); ];
  %wfa = [ q_d * ( 2*rand( n, 1 ) - 1 ); zeros( n, 1 ); ];
  
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
