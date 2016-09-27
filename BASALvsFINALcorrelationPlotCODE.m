
close all
%for timecourse plots
tshift = 8;

cd /Users/frick/Documents/MATLAB/AllImagingDataCompiled
% load('tracesstructdataforPlotting.mat') %2015_11_16
load('ScalarStructure20160229.mat');
load('tracesstructabs-4.mat');%2016_02_19
load('tracesstructfc-4.mat');%2016_02_19


% make a plot of cells that are still tracked after 120 minutes
ds = {'i2dot40i'};
i=1;
tinterval=4;
frame=8;
ds = char(ds);

abstraces = tracesstruct;
abscellies = abstraces.(ds);


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
r = [];
p = [];

        abscellies(abscellies==0) = NaN;
        
        time = 1:1:size(abscellies,2);
        timeadj = (time-frame).*tinterval;

tim = (time-6).*4;
tim = timeadj;
cycle=0;
tpoints = [];
for i = 7:8
    for j = 16:22
        cycle=cycle+1;
        off = abscellies(:,i);
        on = abscellies(:,j);
        indi = ~isnan(off) & ~isnan(on);
        off = off(indi);
        on = on(indi);
        [r(cycle),p(cycle)] = corr(off,on,'type','Pearson','rows','all');
        
        tpoints{cycle} = strcat(num2str(tim(i)),'x',num2str(tim(j)),'-','r=',num2str(round(r(cycle),2,'significant')),'N=',num2str(sum(indi)));
        mdl = fitlm(off,on);
        subplot(4,4,cycle);plot(mdl);title(tpoints{cycle})
        
    end
end



figure(2)
plot(timeadj,abscellies)

figure(3)
%time = [-4,32]
        off = abscellies(:,7);
        on = abscellies(:,16);
        
mdl = fitlm(off,on);
mdl
subplot(2,2,1);plot(mdl); title('lin vs lin')

mdl = fitlm(log(off),log(on));
mdl
subplot(2,2,2);plot(mdl);title('log vs log')


mdl = fitlm(log(off),on);
mdl
subplot(2,2,3);plot(mdl);title('log vs lin')

mdl = fitlm(off,log(on));
mdl
subplot(2,2,4);plot(mdl);title('lin vs log')

        

