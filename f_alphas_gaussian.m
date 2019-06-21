function x = f_alphas_gaussian( n, q_d, alphas,cutfreq,fs)
%   X = F_ALPHAS_GAUSSIAN( N, Q_D, ALPHAS,CUTFREQ,FS)
%
%   Generates a discrete multi-colored (gaussian) noise vector of size N with power 
%   spectrum distribution (in order) of 
%   
%   ALPHA = [a_1 a_2 ... a_M]; 
%
%   all a_i>0, and change of 'color' at frequencies
%
%   CUTFREQ = [f_1 f_2 ... f_(M-1)];
%
%   all must be less than FS. The white noise filtered with these color
%   spetra by bands is Gaussian with (0,Q_D) distribution
%
%   Example: >> alphas = [.8 2/3 1/3]; cutfreq = [10 800];fs=10000;n=1e7;
%            >> q_d = 2;  
%            >> x = f_alphas_gaussian(n,q_d,alphas,cutfreq,fs);
%     

%%  Set the deviation of the noise.
    q_d = sqrt ( q_d );

%%  Generate the coefficients Hk.
    hfa = hfilter(n,alphas);
    wfa = [ q_d * randn( n, 1 ); zeros( n, 1 ); ];
    
    fr = round(cutfreq/fs*n)/n;hfp = [];
    for k=1:length(cutfreq)
        if k==1
            f = [0 fr(1) fr(1) 1]; 
            m = [1 1 0 0];
        else
            f = [0 fr(k-1) fr(k-1) fr(k) fr(k) 1]; 
            m = [0 0 1 1 0 0];   
        end
        hfp = [hfp; fir2(n,f,m,n,10)]; 
    end
    f = [0 fr(end) fr(end) 1]; 
    m = [0 0 1 1 ];
    hfp = [hfp; fir2(n,f,m,n,10)];
    hfp = hfp';hfp = [hfp;zeros(n-1,length(cutfreq)+1)];
%%  Perform the discrete Fourier transforms of Hk and Wk.
    fw= fft(wfa);
    fh= fft(hfa);
    fp= fft(hfp);clear wfa hfa hfp f m q_d fs cutfreq;%also deleting unnecessary data
    corrf=1;
    for k=1:length(fr),
        corrf=abs(fh(fr(k)*n+1,k))/abs(fh(fr(k)*n+1,k+1));
        fh =[fh(:,1:k), corrf*fh(:,k+1:end)];
    end
    fh = fh.*fp;clear fp fr k;
    fh = sum(fh,2);
%%  Multiply the two complex vectors.
    fh = fh(1:n+1,:);fw = fw(1:n+1,:);
    fw = fh .* fw;clear fh;
%  This scaling is introduced only to match the behavior
%  of the Numerical Recipes code...
    fw(1)       = fw(1) / 2;
    fw(end)     = fw(end) / 2;
%%  Take the inverse Fourier transform of the result.
    fw = [ fw; zeros(n-1,1); ];
    x = ifft( fw );
    x = 2*real( x(n/2:3*n/2-1));

  return
  
  function hfa = hfilter(n,alphas)
    dim = length(alphas);
    
    hfa = zeros (2*n, dim );
    hfa(1,:) = ones(1,dim); 
    for i = 2:n
        hfa(i,:) = hfa(i-1,:) .* ( 0.5 * alphas + ( i - 2 ) ) / ( i - 1 );
    end
  end
end