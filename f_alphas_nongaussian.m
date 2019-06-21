function X = f_alphas_nongaussian(N,fs,cdf,cdf_var,colour,colour_var)
%   X = F_ALPHAS_NONGAUSSIAN(N,FS,CDF,CDF_VAR,COLOUR,COLOUR_VAR)
%
%   This function produces a series with length N and sampling rate FS where
%
%   CDF is the probability asociated to the stochastic process and
%
%   CDF_VAR the variables defining the distribution. There are four (4)
%   options:
%
%   1) Normal: cdf = 'Normal';cdf_var = [mean variance];
%
%   2) Uniform: cdf = 'Uniform';cdf_var = [mean variance];
%
%   3) Poisson: cdf = 'Poisson';cdf_var = [mean];
%
%   4) alpha Stable: cdf = 'aStable';cdf_var = [alpha beta mean variance];
%   in this case beta=1 gives Levy, and beta=2 Cauchy (-1<=beta<=1).
%   Besides 0<alpha<=2, and alpha=1/2 is for a Levy distribution.
%
%   COLOUR defines the spectrum associated to the process, it is defined by
%   a string of what you would write in the command window with variable X,
%   e.g.: 'a./X.^b'. The values of the variables are stored in 
%
%   COLOUR_VAR as a vector (in order), e.g., [a b].
%
%
Sxx = spectrumf(colour,colour_var,N,fs);clear colour colour_var;
Xf  = sqrt(2*N*pi*fs*Sxx);
Xf  = ifftshift(Xf);
X   = rnd_cdf(cdf,cdf_var,N); var_Sxx= 2*pi*fs*sum(Sxx)/(N-1); clear cdf cdf_var;
X   = X*sqrt(var_Sxx/var(X)); mean_X =mean(X); clear Sxx fs;
X   = X - mean_X;
[X0, ~ ] = sort(X);
k=1; indX0 = zeros(1,N);
while(k)
    X_phase = angle(fft(X));
    X = real(ifft(exp(1i.*X_phase).*abs(Xf)));
    [~, indX] = sort(X);
    X(indX)=X0;k=k+1;
    if indX==indX0, k=0;end
    indX0 = indX;
end
X = X + mean_X;
end

function rnd = rnd_cdf(cdf,cdf_var,N)
    switch cdf
        case 'Normal'
            rnd = random(cdf,cdf_var(1),cdf_var(2),1,N);
        case 'Uniform'
            rnd = random(cdf,cdf_var(1),cdf_var(2),1,N);
        case 'Poisson'
            rnd = random(cdf,cdf_var,1,N);
        case 'aStable'
            rnd = stblrnd(cdf_var(1),cdf_var(2),cdf_var(4),cdf_var(3),1,N);
     end
end

