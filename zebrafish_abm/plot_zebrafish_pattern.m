%Code adapted from Volkening & Sandstede (2018)
%Plot the simulated zebrafish pattern
% savefolder: folder where simulation outcome is saved
% simname: the name of the simulation (prefix of .mat file)
% t: time point to plot pattern
% savemode:
%  0: display the figure, do not save as png
%  1: display the figure, save as png
%  2: display the figure, save as png and eps
%  3: do not display, do not save (return figure handle only)
%  4: do not display, save the figure as png
%  5: do not display, save the figure as png and eps
%  savemode is an axis handle: plot on that axis, do not save
% crop: 0: show full sized domain. 1: crop to show the most interesting
% part of the pattern, and cut out the figure margins

function [fig] = plot_zebrafish_pattern(savefolder,simname,t,savemode,crop)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
file=strcat(savefolder,'/',simname,'.mat');
if ~exist('crop','var')
    crop = 0; % default setting is no crop
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading simulation results for plotting
load(file,'cellsM', 'cellsXc','cellsXsn','cellsId','cellsIl', 'numMel', 'numXansn', 'numXanc','numIrid','numIril','boundaryY','boundaryX')


%%% Remanming variables
melN = numMel;
xanNc = numXanc;
xanNsn = numXansn;
iriNd = numIrid;
iriNl = numIril;

timeOfInterest=t;

%%% Defining matrices to hold cells on the given day of interest
iriD = cellsId(1:iriNd(timeOfInterest),:,timeOfInterest);
iriL = cellsIl(1:iriNl(timeOfInterest),:,timeOfInterest);
mel = cellsM(1:melN(timeOfInterest),:,timeOfInterest);
xanC = cellsXc(1:xanNc(timeOfInterest),:,timeOfInterest);
xanSN = cellsXsn(1:xanNsn(timeOfInterest),:,timeOfInterest);

if isa(savemode,'matlab.graphics.axis.Axes')
    axes(savemode);
    fig=0;
else
    if savemode>2
        fig = figure(Visible="off");
    else
        fig = figure(Visible="on");
    end
end
hold on

plot(mel(:,1), mel(:,2), 'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize', 3)

plot(iriL(:,1), iriL(:,2), 'v',...
    'MarkerEdgeColor',[58/255 95/255 205/255],...
    'MarkerSize', 4)

plot(iriD(:,1), iriD(:,2), 's',...
    'MarkerEdgeColor',[192/255 192/255 192/255],...
    'MarkerSize', 4)

plot(xanC(:,1), xanC(:,2), '*',...
    'MarkerEdgeColor',[255/255  185/255 15/255],...
    'MarkerFaceColor',[255/255  185/255 15/255],...
    'MarkerSize', 3)

plot(xanSN(:,1), xanSN(:,2), 'yx','MarkerSize',2)

if ~crop
    %%% scale-bar
    plot([boundaryX(timeOfInterest)-600 boundaryX(timeOfInterest)-100],[100 100],'r-','linewidth',5, 'color',[205/255 0/255 0/255]);

    % box around pattern
    plot(-50:1:boundaryX(timeOfInterest)+50, ones(1,boundaryX(timeOfInterest) +101)*(boundaryY(timeOfInterest)+50),'k','LineWidth',1.5)
    plot(-50:1:boundaryX(timeOfInterest)+50, -50*ones(1,boundaryX(timeOfInterest)+101),'k','LineWidth',1.5)
    plot(-50*ones(1,boundaryY(timeOfInterest)+101), -50:1:boundaryY(timeOfInterest)+50,'k','LineWidth',1.5)
    plot((boundaryX(timeOfInterest)+50)*ones(1,boundaryY(timeOfInterest)+101), -50:1:boundaryY(timeOfInterest)+50,'k','LineWidth',1.5)

    k1=strfind(file,'/');
    k1=k1(end);
    title(sprintf("Sim id %s, Day %d (%d dpf)",file(k1+1:end),timeOfInterest,timeOfInterest+20),'FontSize',14,'Interpreter', 'none');
    xlabel('fish length x');
    ylabel('fish width y');
    xlim([0 4000]);
    ylim([-1000,3000]);
else
    %%% scale-bar
    plot([boundaryX(timeOfInterest)-600 boundaryX(timeOfInterest)-100],[150 150],'r-','linewidth',5, 'color',[205/255 0/255 0/255]);

    fig.Position = [0,0,500,500/boundaryX(timeOfInterest)*(boundaryY(timeOfInterest)-100)];
    axis equal;
    xlim([0,boundaryX(timeOfInterest)]);
    ylim([50,boundaryY(timeOfInterest)-50]);
    xticks([]);
    yticks([]);
    tightEdge(fig);
    %pbaspect([1 1 1]);
end

hold off

if ~isa(savemode,'matlab.graphics.axis.Axes')
    if (savemode==1) || (savemode==2) || (savemode==4) || (savemode==5)
        saveas(fig,sprintf('%s/%s_pattern_t=%04.1f.png',savefolder,simname,t),'png');
    end
    if (savemode==2) || (savemode==5)
        saveas(fig,sprintf('%s/%s_pattern_t=%04.1f.eps',savefolder,simname,t),'epsc');
        saveas(fig,sprintf('%s/%s_pattern_t=%04.1f.svg',savefolder,simname,t));
        saveas(fig,sprintf('%s/%s_pattern_t=%04.1f.fig',savefolder,simname,t));
    end
end

end