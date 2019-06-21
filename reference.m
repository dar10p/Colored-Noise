n = 100000; alpha1 =2/3;q_d=1;alpha2 =1/3;

hfa1 = zeros (2*n, 1 );
  hfa1(1) = 1.0; 
  for i = 2:n
    hfa1(i) = hfa1(i-1) * ( 0.5 * alpha1 + ( i - 2 ) ) / ( i - 1 );
  end
    hfa1(n+1:2*n) = 0.0;
 
hfa2 = zeros (2*n, 1 );
  hfa2(1) = 1.0; 
  for i = 2:n
    hfa2(i) = hfa2(i-1) * ( 0.5 * alpha2 + ( i - 2 ) ) / ( i - 1 );
  end
    hfa2(n+1:2*n) = 0.0;
    wfa = [ q_d * randn( n, 1 ); zeros( n, 1 ); ];
%%
fs =10; cutfreq=2.5;fr = round(cutfreq/fs*n)/n;
    f = [0 fr fr 1]; m = [1 1 0 0];
    hfb = fir2(n,f,m,n,10); hfb = hfb';
    m = [0 0 1 1];
    hfc = fir2(n,f,m,n,10); hfc = hfc';
%%
%     [fhb,w]=freqz(hfb,1,2*n,'whole');
%     [fhc,w]=freqz(hfc,1,2*n,'whole');  
%     [fha1,w]=freqz(hfa1,1,2*n,'whole');
%     [fha2,w]=freqz(hfa2,1,2*n,'whole');
%     [fw,w]=freqz(wfa,1,2*n,'whole');
%%
hfb = [hfb; zeros(n-1,1)];
hfc = [hfc; zeros(n-1,1)];
%%
ifhb = fft(hfb);
ifhc = fft(hfc);
ifha2 = fft(hfa2);
ifha1 = fft(hfa1);
ifw = fft(wfa);
%%
corrf=abs(ifha1(fr*n+1))/abs(ifha2(fr*n+1));
  ifh1 = ifhb.*ifha1;ifh2=corrf*ifhc.*ifha2;
%%
%   corrf=abs(fha1(fr*n+1))/abs(fha2(fr*n+1));
%   fh1 = fhb.*fha1;fh2=corrf*fhc.*fha2;
%%
 ifh1 = ifh1(1:n+1);ifh2 = ifh2(1:n+1);
 ifw = ifw(1:n+1);
      %%
 ifw = (ifh1+ifh2).*ifw;

%%

    ifw(1)       = ifw(1) / 2;
    ifw(end)     = ifw(end) / 2;
    ifw = [ifw;zeros(n-1,1); ];
  %%
  x = ifft( ifw );
  
  x = 2*real( x(n/2:3*n/2) );