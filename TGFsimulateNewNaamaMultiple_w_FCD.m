function TGFsimulateNewNaamaMultiple_w_FCD(X,protein)
global Tgfz Tgfon Tgfoff Tgfbasal 
close all

%initialize structures and figures
BB = struct();
CC = struct();
DD = struct();
firstfigure=4;
secondfigure=5;
%load conditions
c = feval('TGFconditions_w_FCD');     
save conditions_beta.dat c -ascii;
%load parameters
p = feval('TGFparameters_w_FCD');
save parameters_beta.dat p -ascii;

%Dynamics values
tn    = c(1);     % Time span for integration,seconds
Tgfoff = c(2);
Tgfbasal = c(3);
Tgfon = c(4);
% tspan = [0:5:tn];
tspan = [0;tn];
Tgfz = 's';

number_of_doses=2;
if number_of_doses == 1
    TgfF = 1;
else
TgfF = log10(logspace(0.02,1,number_of_doses));
end

for FFF = 1:length(TgfF)-1
% for FFF = 5
Tgfon = TgfF(FFF+1);
%==========================================================================
% Computing the unperturbed and perturbed solutions (dimensional solution)
%==========================================================================
% Increasing one parameter at a time by X-fold
%--------------------------------------------------------------------------
setvar = 0.1;
COLORS = colormap(parula(11));
iterations = 1;
% PertToEmploy = [1 ; X ; 1./X]; %[unperturbed; perturbed; perturbed down];
for i=1:iterations;        % 11 parameters
    disp(i) 
%make parameter perturbation
    p = feval('TGFparameters_w_FCD');
    variation = lognrnd(0,0.1,length(p),1);
    variation(9) = lognrnd(0,0.4,1,1);
    p = p.*variation';
    save parameters_beta.dat p -ascii;

%============
%Time course for basal state
%============
    y0 = feval('TGFconcentrations_w_FCD');
    save concentrations.dat y0 -ascii;
    y0(22) = Tgfoff;
    %Computing initial guess for dy, using decic
    fixed_y0 = ones(size(y0));
    fixed_dy0 = zeros(size(y0));
    dy0 = zeros(size(y0));
    [y0mod,dy0mod] = decic('TGFequations_w_FCD',0,y0,fixed_y0,dy0,fixed_dy0);
  % Solving the ODEs
    [TT,YY] = ode15i('TGFequations_w_FCD',tspan,y0mod,dy0mod);

%============
%Time course for stimulated state
%============    
    y0 = YY(end,:);    
    y0(22) = Tgfon;
 %  Computing initial guess for dy, using decic
    fixed_y0 = ones(size(y0));
    fixed_dy0 = zeros(size(y0));
    dy0 = zeros(size(y0));
    [y0mod,dy0mod] = decic('TGFequations_w_FCD',0,y0,fixed_y0,dy0,fixed_dy0);
  % Solving the ODEs
    [TTT,YYY] = ode15i('TGFequations_w_FCD',tspan,y0mod,dy0mod);

%============
%Concatentate the time courses for both states
%============    
T = vertcat(TT,TTT+TT(end));
Y = vertcat(YY,YYY);

figure(2222)
subplot(2,1,1);
plot(T,Y(:,3));hold on
plot(T,Y(:,24));hold on
plot(T,Y(:,25));hold on


alp = find(Y(:,25) == max(Y(:,25)),1,'last');
alpi = find(T < 4.79e05,1,'last');
alpt = find(T > 4.95e05,1,'first');

initi = alpi;
termini = alpt;

subplot(2,1,2);
plot((T(initi:termini)-4.8e05)/60,Y(initi:termini,3)./((Y(initi,3))));hold on
plot((T(initi:termini)-4.8e05)/60,Y(initi:termini,24)./((Y(initi,24))));hold on
plot((T(initi:termini)-4.8e05)/60,Y(initi:termini,25)./((Y(initi,25))));hold off
xlim(([4.79e05 5.05e05]-4.8e05)./60)


totalSmad = p(9);
BB(i).Expression = totalSmad;
BB(i).TimeCourse = Y;
BB(i).Time = T;
% BB(i).PerturbationStrength = prod(variation);
BB(i).PerturbationStrength = 10.^sum(abs(log10(variation)));
BB(i).Parameters = p;
BB(i).Dose = Tgfon;
% BB(i).Color = COLORS{perturbedParameter};


%% species  
details = chooseSpecies(T,Y);
species = details.(protein);
totalspecies = details.S2total;
basal = find(T<max(tspan),1,'last');

FCT = zeros(size(Y));
for j = 1:size(Y,2);
FCT(:,j) =  Y(:,j)./Y(basal,j);
BB(i).FCTimeCourse = FCT;
end

% % if PP ==1
% figure(firstfigure)
% h = subplot(2,1,1);plot((T-tn)./60,species);hold on
% h.XTick = -120:60:600;
% h.XTickLabel = -120:60:600;
% xlim([-120 600])
% h =subplot(2,1,2);plot((T-tn)./60,species./(species(basal)));hold on
% h.XTick = -120:60:600;
% h.XTickLabel = -120:60:600;
% xlim([-120 600])
% % end
%% descriptors
CC = datastructmaker(protein,T,Y,basal,totalspecies,i,CC);
DD = datastructmaker('S24nuc',T,Y,basal,totalspecies,i,DD);
end
%% figures    
% figure(secondfigure)
% scatter(vertcat(BB.PerturbationStrength),vertcat(CC.peak)./nanmedian(vertcat(CC.peak)));hold on
% scatter(vertcat(BB.PerturbationStrength),vertcat(CC.foldchange)./nanmedian(vertcat(CC.foldchange)));hold on
% scatter(vertcat(BB.PerturbationStrength),vertcat(CC.percen)./nanmedian(vertcat(CC.percen)));hold on
% scatter3(vertcat(BB.Dose),vertcat(DD.percen),vertcat(CC.percen));hold on
BBB{FFF}=BB;
CCC{FFF}=CC;
DDD{FFF}=DD;


fnames = fieldnames(CC)';
for jjj=1:length(fnames)
InfoMAT(jjj,:) = vertcat((CC.(fnames{jjj})));
end
INFOyo{FFF} = InfoMAT;
end
% plotSensitivityAnalysis(BB,p)
stophere=1;
save('/Users/frick/Documents/Goentoro_Lab/DATA/Modeling/2015_08_24 Supplement Modeling/Naama p vary 01 while S3 vary 04/CONDITIONS/DataAfterTGFsimulateNEW.mat')


% for i=1:length(TgfF)
% IMat = INFOyo{i};
% for j = 1:size(IMat,1)
% subplot(3,3,j);scatter(randi(1000,[1 1000]),IMat(j,:));hold on
% end
% end

end


function quant = quantifyThem(Time,Basal,species)
quant = 1;
end


  function details = chooseSpecies(T,Y)
details.S2nuc = Y(:,11)+Y(:,13)+Y(:,16)+Y(:,18); %nuclear Smad2
details.S2cyto =  Y(:,1)+Y(:,3)+Y(:,6)+Y(:,8); %cytoplasmic Smad2
details.S2total = Y(:,11)+Y(:,13)+Y(:,16)+Y(:,18)+Y(:,1)+Y(:,3)+Y(:,6)+Y(:,8); %total smad2

details.S4nuc = Y(:,15)+Y(:,16); %nuclear S4
details.S4cyto = Y(:,5)+Y(:,6); %cytoplasmic S4
details.S4total = Y(:,15)+Y(:,16)+Y(:,5)+Y(:,6); %total S4

details.S24nuc = Y(:,16); %S24 nuclear
details.S24cyto = Y(:,6);
details.S24total = Y(:,16)+Y(:,6);
  end
    
  

  
  function plotSensitivityAnalysis(BB,p)
  tights = {'maxrate','foldchange','basal','peak','relative rate','percen','NT'};
 xticklabel = {'kExport','kImport','kPhos','kOn','kOff','CIF','kDephos','[PPase]','[Smad3]t','[Smad4]t','kTgfB','',};
 xtick = 1:(length(p)+1).*2;
 xticklabels = horzcat(xticklabel,xticklabel);
 
 
fnames = fieldnames(BB);
    for bb = 1:length(fnames)
    PLOTTER = BB.(fnames{bb});
    
    
    upp = PLOTTER(2,:);
    downn = PLOTTER(3,:);
    
 
    subplot(4,3,bb);h = bar([1:length(upp)],upp,'FaceColor',[0.8 0.85 0.9],'EdgeColor','none');title(fnames{bb});hold on
    subplot(4,3,bb);hh = bar([1:length(upp)]+1+length(ones(size(upp))),downn,'FaceColor',[0.2 0.3 0.2],'EdgeColor','none');title(fnames{bb});hold off
    xlim([0 24])
    
    
    
    m = gca;
    m.YScale = 'log';
    set(m,'XTickLabel',xticklabels,'XTick',xtick);
    h.BaseValue = PLOTTER(1,1);
    hh.BaseValue = PLOTTER(1,1);
    
    ax = h.Parent;
    ax.XTickLabelRotation = 45;
    prim = logspace(log10(PLOTTER(1,1)./5),log10(PLOTTER(1,1).*5),6);
    sec = round(prim,2,'significant');
    tert = sort(horzcat(sec,PLOTTER(1,1)));
    ax.YTick = tert;
    ax.YTickLabel = tert;
    ylim([min(tert) max(tert)]);
    
    end
    stophere=1;
    
  end
    
  
function CC = datastructmaker(protein,T,Y,basal,totalspecies,i,CC)
details = chooseSpecies(T,Y);
species = details.(protein);
rate = gradient(species);% rate
CC(i).maxrate = max(rate(basal+1:end));%max rate
CC(i).maxrelrate = max(rate(basal+1:end))./species(basal);% relative rate
CC(i).foldchange = species(length(species))./species(basal);% fold change
CC(i).basilico = species(basal);% basal
CC(i).peak = species(length(species));% peak
CC(i).percen = CC(i).foldchange-1;%percent
CC(i).NT = species(length(species))./totalspecies(length(totalspecies));% nuclear/total
end

    
