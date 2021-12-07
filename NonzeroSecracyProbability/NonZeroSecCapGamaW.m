clc
clear
close all
%~~~~~~~~~~~~~~~~~~~Non Zero Secrecy Capacity~~~~~~~~~~~~~~~~~~%
% In this code, non zero secrecy capasity in a relay system is simulated.  
%%
N              = 5; %number of relays.
ittr           = 1000000; %number of itterations.
sig            = 1;  %noise sigma
gama_SEdB      = [-5,0,5]; %dB % Source to eavesdropper channel gains
gama_SDdB      = -10:10;    %dB % Source to destination channel gains
gama_SD        = db2pow(gama_SDdB); %Transform from dB to power.
marks          = ['s';'*';'o';'d';'>']; %Each marks is used for one plot.
color_codes    = [0,0,0.9;0.8,0,0.8;0.95,0,0;]; % Codes for diffetent color of plots.
Rho            = 0.5;  %Delay parameter 
rho            = [1 sqrt(Rho);sqrt(Rho) 1];
chol_s         = chol(rho); %This parameter is used to model the correlated channel in outdated relay selection.
%%
% In this loop, non zero secrecy capacity is calculated and ploted...
% for different amount of gamma_SEdB.
for k=1:length(gama_SEdB)
    gama_se             = db2pow(gama_SEdB(k));
    SecrecyCapacityTh   = zeros(length(gama_SD),1);
    SecrecyCapacityMont = zeros(length(gama_SD),1);
    for j=1:length(gama_SD)
        %%
        SecrecyCapacityTh(j)= NZSCTH(N,gama_SD(j),Rho,gama_se); %This function calculates the non zero secrecy capacity using the the derived formula.
        %%  
        Hsd             = (1/sqrt(2)).*(normrnd(0,sig,1,ittr) + 1i*normrnd(0,sig,1,ittr)); %source ro destination channels vector
        Hsr             = (1/sqrt(2)).*(normrnd(0,sig,N,ittr) + 1i*normrnd(0,sig,N,ittr)); % source to relays channels matrix.
        Hrd             = (1/sqrt(2)).*(normrnd(0,sig,N,ittr) + 1i*normrnd(0,sig,N,ittr)); % relays to source channels matrix.
        %Opportunistic relay selection:
        gama_sr         = ((abs(Hsr)).^2)*gama_SD(j); %Source to relays gamma
        gama_rd         = ((abs(Hrd)).^2)*gama_SD(j); %Relays to destination gamma
        min_gama        = min(gama_sr,gama_rd);
        max_min_gama    = max(min_gama);           %Opportunistic relay selection.
        % In this loop best relays are chosen.
        for i=1:ittr
            Rb          = find(min_gama(:,i)==max_min_gama(i)); %Choosing best Relays.
            Hsrc(1,i)   = Hsr(Rb,i); %channel vector of source to best relay.
            Hrdc(1,i)   = Hrd(Rb,i);  %channel vector of best relay to destination. 
        end
        %Chanel coefficient to be used for correlated delayed channels.
        Hsrcc           = (1/sqrt(2)).*(normrnd(0,sig,1,ittr) + 1i*normrnd(0,sig,1,ittr));
        Hrdcc           = (1/sqrt(2)).*(normrnd(0,sig,1,ittr) + 1i*normrnd(0,sig,1,ittr));
        % channel coefficients of compensated delayed channels.
        Hsrto           = chol_s'*[Hsrc;Hsrcc];
        Hrdto           = chol_s'*[Hrdc;Hrdcc]; % Channel at the time of transmition
  
        gama_SRto       = ((abs(Hsrto(2,:))).^2)*gama_SD(j);
        gama_RDto       = ((abs(Hrdto(2,:))).^2)*gama_SD(j);
        gama_Relay      = min(gama_SRto,gama_RDto);
        gamasd_direct   = ((abs(Hsd)).^2)*gama_SD(j);
        GamaSD_updated = gamasd_direct+gama_Relay;   %%%%%%%%%%%%%%%%%%% Updated source to destination gama
        %%
        % Source to eavesdroper and relays to eavesdroper channels coefficients and gamas
        
        Hse             = (1/sqrt(2)).*(normrnd(0,sig,1,ittr) + 1i*normrnd(0,sig,1,ittr));
        gama_se1        = gama_se.*((abs(Hse)).^2);
        HRe             = (1/sqrt(2)).*(normrnd(0,sig,1,ittr) + 1i*normrnd(0,sig,1,ittr));
        gama_RE         = ((abs(HRe)).^2)*gama_se;
        HsR             = (1/sqrt(2)).*(normrnd(0,sig,1,ittr) + 1i*normrnd(0,sig,1,ittr));
        gama_SR         = ((abs(HsR)).^2)*gama_se;
        gama_RelayEve   = min(gama_SR,gama_RE);
        Gama_SE         = gama_se1+gama_RelayEve; % Source to eavesdroper gama.
        %% Computing Nonzero secrecy capacity 
        SecrecyCapacity = 0.5*(log2(1+GamaSD_updated)-log2(1+Gama_SE));
        Num_NzSC_Transmits= length(find(SecrecyCapacity>0));
        SecrecyCapacityMont(j)         = Num_NzSC_Transmits/ittr;
    end
    %% Plot 
    semilogy(gama_SDdB,SecrecyCapacityTh,'color',color_codes(k,:),'LineWidth',1);
    hold on
    semilogy(gama_SDdB,SecrecyCapacityMont,marks(k),'color',color_codes(k,:));
    xlabel('\gamma_{SD}')
    ylabel('Nonzero Secrecy Probability')
end
legend('theory,\gamma_{SE} = -5','simulation,\gamma_{SE} = -5','theory,\gamma_{SE} = 0','simulation,\gamma_{SE} = 0','theory,\gamma_{SE} = 5','simulation,\gamma_{SE} = 5','Location','southeast')
grid