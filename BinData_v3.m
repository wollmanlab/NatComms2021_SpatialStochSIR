function [Xb, Yb, stdXb, stdYb, steXb, steYb,Xts, Yts] = BinData_v3(X, Y, bins, varargin);

nanzero = ParseInputs('nanzero',false,varargin);

if nanzero
    Y(Y==0)=nan;
end

if numel(bins)==1;
    binedges = quantile(X,0:(1/bins):1);
else
    binedges = bins;
end


Xts=cell(1,numel(binedges)-1);
Yts=cell(1,numel(binedges)-1);
for i=1:(numel(binedges) -1),
    Xt = X(logical((X<binedges(i+1)) .* (X>=binedges(i))));
   
       Xts{i} = Xt;
    Xb(i) = nanmean(Xt);
    stdXb(i) = nanstd(Xt);
    steXb(i) = stdXb(i)/sqrt(length(Xt));
    
    Yt = Y(logical((X<binedges(i+1)) .* (X>=binedges(i))));
    Yts{i} = Yt;
    Yb(i) = nanmean(Yt);
    stdYb(i) = nanstd(Yt);
    steYb(i) = stdYb(i)/sqrt(length(Yt));
end;