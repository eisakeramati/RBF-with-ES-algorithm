function [ind,maxi]=indmax(x)
% [ind,maxi]=indmax(x)
% Finds index to and maximum of N-dimensional array
% ind is cell array
%
% Example:
% r=rand([2 3 4])
% [ind,maxi]=indmax(r)
% so that
% maximum=r(ind{:})
[x1,i1]=max(x);
x1=squeeze(x1);
i1=squeeze(i1);
if numel(x1)>1
    [i2,x2]=indmax(x1);
    ind={i1(i2{:}),i2{:}};
    maxi=x2;
else
    ind={i1};
    maxi=x1;
end