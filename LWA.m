% function [Y,UMFYy,UMFYmu,LMFYy,LMFYmu] = LWA(X,W,n)

%
% [Y, UMFYy,UMFYmu,LMFYy,LMFYmu] = LWA(X,W,n)
%
% To compute the LWA [1,2] for IT2 FSs determined by the nine 
% parameters in Fig. 1 of Readme.doc.
%
% [1] Dongrui Wu and Jerry M. Mendel, “Aggregation using the linguistic 
% weighted average and interval type-2 fuzzy sets,” IEEE Trans. on Fuzzy 
% Systems, vol. 15, no. 6, pp. 1145--1161, 2007.
%
% [2] Dongrui Wu and Jerry M. Mendel, “Corrections to ‘Aggregation using 
% the linguistic weighted average and interval type-2 fuzzy sets’,” IEEE 
% Trans. on Fuzzy Systems, in press.
%
% Dongrui WU (dongruiw@usc.edu), 11/16/2008
%
% X and W: MFs of the subcriteria and weights. They have the same number of
% rows.
% n: number of alpha-cuts. Default is 2.
%
% Y: the LWA approximated by 9 parameters.
% UMFYy and UMFYmu: y- and mu-coordinates of the UMF of the LWA
% LMFYy and LMFYmu: y- and mu-coordinates of the LMF of the LWA

X = [0 0 1 3 0 0 1 2 0.8; 0 2 3 5 1 2 3 4 0.8; 2 4 5 5 3 4 5 5 0.8];

W = [0 0 0 0.4 0 0 0 0.2 1; 0.2 0.6 0.6 0.8 0.4 0.6 0.6 0.6 1; 0.6 1 1 1 0.8 1 1 1 1];

n = 2;

% if nargin==2 %% set default n
%     n=2;
% end

if size(X,2)==8
    X(:,9)=1; % set default height
end
if size(W,2)==8
    W(:,9)=1; % set default height
end

[Yu,UMFYy,UMFYmu] = FWA(X(:,1:4),W(:,1:4),n);
[Yl,LMFYy,LMFYmu] = FWA(X(:,5:9),W(:,5:9),n);
Y=[Yu(1:4) Yl];

plotIT2(UMFYy,UMFYmu,LMFYy,LMFYmu);
