function cmap = redblue(n)
% redblue  Diverging red-white-blue colormap.
if nargin < 1, n = 256; end
half = floor(n/2);
bot  = [linspace(0.12,1,half)', linspace(0.53,1,half)', linspace(0.90,1,half)'];
top  = [linspace(1,0.70,n-half)', linspace(1,0.09,n-half)', linspace(1,0.09,n-half)'];
cmap = [bot; top];
end
