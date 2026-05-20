%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting utility for persistent homology barcodes.
%
%
%Copyright (C) Yue Liu, Alexandria Volkening, 2026
%
%This program is free software: you can redistribute it and/or 
%modify it under the terms of the GNU General Public License as 
%published by the Free Software Foundation, version 3.
%
%This program is distributed in the hope that it will be useful, 
%but WITHOUT ANY WARRANTY; without even the implied warranty of 
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public 
%License along with this program. If not, see 
%https://www.gnu.org/licenses/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

