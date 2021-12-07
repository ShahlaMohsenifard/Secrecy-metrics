% This function calculates the formula of outage probability which is 
%derived in the paper.  
function  [P_O_Th] = Out_Prob_Th(N,gama_SD,Rho,gama_se,Rs)
gama_ce         = gama_se/2;
gama_c          = gama_SD/2;
P_out           = zeros(N,1);
Rs=2*Rs;
for n=1:N
    b           = (gama_c*(n*(1-Rho)+Rho))/(1.5*n);
    nN          = nchoosek(N,n);
    P_out(n)    = nN*((-1)^(n-1))*(1/(gama_ce-gama_se))*(1/(b-gama_SD))*...
        (b*exp((1-2^Rs)/b)*((1/(2^Rs/b+1/gama_ce))-(1/(2^Rs/b+1/gama_se)))...
        -gama_SD*exp((1-2^Rs)/gama_SD)*((1/(2^Rs/gama_SD+1/gama_ce))-(1/(2^Rs/gama_SD+1/gama_se))));
end
P_O_Th          = 1-sum(P_out);