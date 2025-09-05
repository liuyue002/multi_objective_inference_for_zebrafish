
%Utility for cutting out the margins of a figure
function tightEdge(fig)
axs=fig.Children;
ax=axs(1);
outer = ax.OuterPosition;
ti = ax.TightInset;
l = outer(1)+ti(1);
b = outer(2)+ti(2);
width = outer(3) - ti(1) - ti(3);
height = outer(4) - ti(2) - ti(4);
ax.Position = [l b width height];
end