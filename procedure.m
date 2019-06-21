%% Fit Results
syms z;
P06Bx1 = solve([num2str(PD06B_fittedmodel2.p1) '*z +' num2str(PD06B_fittedmodel2.p2)...
    '=' num2str(PD06B_fittedmodel1.p1) '*z +' num2str(PD06B_fittedmodel1.p2)],z);

%% Plot Results
p06ax3 = x(x>=4.6);
%p06ax2 = p06ax2(p06ax2>=3);
p06ay3 = PD06A_fittedmodel3.p1*p06ax3+PD06A_fittedmodel3.p2;
