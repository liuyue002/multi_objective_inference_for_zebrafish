function [fig] = plot_pland(A,visible)
%Plot persistence landscape. Only the first 5 p-landscape curves are plotted.
% A: a matrix encoding p-landscape curves
% visible: whether to make the plot window visible

fig=figure(Visible=visible);
num_plot = 5;
hold on;
for i=1:min(num_plot,size(A,1))
    plot(A(i,:));
end
hold off;
legend(arrayfun(@line_label,1:num_plot,UniformOutput=0),Interpreter='latex');
xlabel('$m (\mu m)$',Interpreter='latex');
ylabel('$h (\mu m)$',Interpreter='latex');
title('persistence landscape');
end

%%
function str = line_label(n)
str = sprintf('$\\lambda_%d$',n);
end