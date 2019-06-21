function C = spectrumf(colour, colour_var,N,fs)
%   Sxx = SPECTRUMF(COLOUR,COLOUR_VAR,N,FS)
%
%   COLOUR is the continuous spectrum, expresed as a string function of X;
%   that is, 'a./X.^b'. The constant variables must be from the dictionary 
%   'abcdefghijklmnopqrstuvwxyz'.
%
%   COLOUR_VAR are the values given to the constant variables: [.5 5]
%   asigns a =.5; b=5;
%
%   N is the sample size
%
%   FS is the sampling frecuency
%
%
    num_var = length(colour_var);
    dict = 'abcdefghijklmnopqrstuvwxyz';%dictionary for variables
    v = cell(1,num_var);
    for i=1:num_var;
        v{i} = dict(i);
    end
    vars = genvarname(v);%creates variables
    for i=1:num_var
        eval([vars{i} '=colour_var(' num2str(i) ');']);%assign values to variables
    end
    X = (2*pi*fs/N)*((0:N-1)-N/2);%frequencies where the spectrum is evaluated
    C = eval(colour); C(N/2+1)=0;%the spectrum
end