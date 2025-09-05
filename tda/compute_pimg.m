function [fig,X,Y,Z] = compute_pimg(A, tda_params, display_plot, savefile)
%Produce TDA persistence diagram/image
% A is either a N*2 matrix, each row is a barcode with birth/death time
% or a path to a text file encoding such a matrix
% tda_params: a struct that contains the parameters related to TDA, should
% contain following fields:
%  bmax, pmax: max birth/persistence for the diagram
%  num_pts: resolution of the diagram
%  phi: a smoothing kernel (usually a Gaussian)
%   it should have the signature phi = @(x,y,u,v)...,
%   where u,v are the centre of the kernel
%  w: a weighting function w=@(b,p), where b=birth, p=persistence
%   good choices are 1, p, or sigmoid:
%   1/(1+exp(-(p-p0)/pw)), with a inflection pt at p0 and width pw
%   if phi or w are not function handle, produce a plain persistence diagram
%  phimask: a convolution mask computed from phi. Optional.
% plot: whether to display the diagram
% savefile: path for saving the diagram (only if it's a string)
% return handle to the figure, and persistence surface

if (isstring(A)||ischar(A)) && exist(A,'file')
    % A is a file
    A = readmatrix(A,FileType="text");
end
assert(ismatrix(A) && size(A,2)==2, 'invalid input');
num_bar=size(A,1);
A(:,2) = A(:,2)-A(:,1); % transform to birth-persistence matrix

if display_plot
    fig = figure(Visible="on");
else
    fig = figure(Visible="off");
end

bmax = tda_params.pimg_bmax;
pmax = tda_params.pimg_pmax;
num_pts = tda_params.pimg_numpts;
if isfield(tda_params,'phi')
    phi = tda_params.phi;
else
    phi=nan;
end
if isfield(tda_params,'w')
    w = tda_params.w;
else
    w=nan;
end
if isfield(tda_params,'phimask')
    phimask = tda_params.phimask;
else
    phimask = nan;
end

if isa(phi,'function_handle') && isa(w,'function_handle')
    % plot persistence image
    xx=linspace(0,bmax,num_pts);
    yy=linspace(0,pmax,num_pts);
    db=xx(2)-xx(1);
    dp=yy(2)-yy(1);
    [X,Y] = meshgrid(xx,yy);
    Z = zeros(size(X));

    if isnan(phimask)
        xx2=linspace(-bmax,bmax,num_pts*2-1);
        yy2=linspace(-pmax,pmax,num_pts*2-1);
        [X2,Y2] = meshgrid(xx2,yy2);
        phi0=@(x,y) phi(x,y,0,0);
        phimask = arrayfun(phi0,X2,Y2);
    end

    for i=1:num_bar
        b=A(i,1);
        p=A(i,2);
        if (b > bmax) || (p>pmax)
            continue;
        end
        b_ind = round(b/db)+1;
        p_ind = round(p/dp)+1;
        Z(b_ind,p_ind) = Z(b_ind,p_ind) + w(b,p);
    end
    Z = Z';
    Z = conv2(Z,phimask,'same');

    % asssume values smaller than 1e-5*Zmax don't matter
    % the sparsity makes save file a lot smaller
    Zmax = max(Z,[],'all');
    Z(Z<Zmax*1e-5) = 0;


    surf(X,Y,Z,EdgeColor="none");
    xlabel('Birth $(\mu m)$',Interpreter='latex');
    ylabel('Persistence $(\mu m)$',Interpreter='latex');
    if max(Z,[],'all') > 0
        clim([0,max(Z,[],'all')]);
    else
        clim([0,1]);
    end
    view(2);
    xlim([0,bmax]);
    ylim([0,pmax]);
    colorbar;
    %title('Persistence image');
else
    X=nan;
    Y=nan;
    Z=nan;
    hold on;
    for i=1:num_bar
        plot(A(i,1),A(i,2),'*k');
    end
    hold off;
    xlabel('Birth $(\mu m)$',Interpreter='latex');
    ylabel('Persistence $(\mu m)$',Interpreter='latex');
    xlim([0,bmax]);
    ylim([0,pmax]);
    %title('Persistence diagram');
end
if (isstring(savefile)||ischar(savefile))
    saveas(fig,savefile,'png');
end

end