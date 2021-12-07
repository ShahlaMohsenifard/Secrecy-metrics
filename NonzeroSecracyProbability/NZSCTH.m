
%This function calculates the formula of "non zero secrecy capacity", 
%which is derived in the paper. 
function [NZSecCapTh]=NZSCTH(N,gama_SD,Rho,gama_Se)
gama_C         = gama_SD/2;
gama_Ce        = gama_Se/2;
pcs_n          = zeros(N,1);
for n=1:N
    nN          = nchoosek(N,n);
    b           = (gama_C*(n*(1-Rho)+Rho))/n;
    pcs_n(n)    = nN*((-1)^(n-1))*(1-(1/(b-gama_SD))*(1/(gama_Ce-gama_Se))...
        *(gama_Ce*(1/(1/b+1/gama_Ce)-1/(1/gama_SD+1/gama_Ce))-gama_Se*...
        ((1/(1/b+1/gama_Se)-1/(1/gama_SD+1/gama_Se)))));
end
NZSecCapTh=sum(pcs_n);
