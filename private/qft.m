function  Q = qft(n)
%QFT Create Quantum Fourier Transform matrix
%   Q=QFT(N) Create N-bit QFT matrix with M=2^N rows and cols

%  QFT(N)*A = sqrt(length(A))*ifft(A) 
%  QFT(N)'*A = fft(A)./sqrt(length(A)) 

m = 2^n;
w=exp(2*pi*1i/m);
row=0:m-1;
% Q = zeros(m);
% for r=1:m
%     Q(r,:)=row*(r-1);
% end
% Q = w.^Q
Q = fliplr(vander(w.^row));
Q = Q/sqrt(m);