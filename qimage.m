function qimage(m)
%QIMAGE plot with black=-1, gray=0, while=+1
%  QIMAGE(M)

image(64*rescale(m, -1, 1, 0, 1))
colormap gray
end

function m = rescale(m, oldmin, oldmax, newmin, newmax)
%m = rescale(m, oldmin, oldmax, newmin, newmax)
%
% Rescales values in matrix m from old scale to new scale.
% (useful for images)

m = (m-oldmin)*((newmax-newmin)/(oldmax-oldmin)) + newmin;
end