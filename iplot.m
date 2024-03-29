function iplot(X)
%IPLOT 3D plot for complex timeseries, eg. DFT
%   IPLOT(X)


N=length(X);
hold on
for i=1:N
    plot3([i i], [0 imag(X(i))], [0 real(X(i))], 'r');
end
plot3(1:N, imag(X), real(X), 'm');
xlabel('t');
ylabel('im');
zlabel('re');

%axis
ax = zeros(N);
plot3(1:N, imag(ax), real(ax), 'k');
hold off