function [fig] = plot_barcode(A,visible)
%Produce TDA barcode plot
% A is either a N*2 matrix, each row is a barcode with birth/death time
% or a path to a text file encoding such a matrix

if (isstring(A)||ischar(A)) && exist(A,'file')
    % A is a file
    A = readmatrix(A,FileType="text");
end
assert(ismatrix(A) && size(A,2)==2, 'invalid input');

num_bar=size(A,1);
fig=figure(Visible=visible);
hold on;
for i=1:num_bar
    plot(A(i,:),[i,i],'-k');
end
hold off;
xlabel('$\epsilon (\mu m)$',Interpreter="latex");
end

