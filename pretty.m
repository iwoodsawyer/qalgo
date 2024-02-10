function str = pretty(struct, b_dec)
%PRETTY Get pretty-printed string depicting a state struct
%   STR=PRETTY(STRUCT) print state struct using binary representation.
%
%   STR=PRETTY(STRUCT,1)  print state struct using decimal representation.

if nargin<2,b_dec=0; end

if ~isstruct(struct) %if vec
    struct(abs(struct)<sqrt(eps))=0;
    struct=vec2struct(struct);
end

str=[];

for i=1:length(struct)
    if (abs(real(struct(i).alpha)) < sqrt(eps))
        struct(i).alpha = imag(struct(i).alpha)*1i;
    end
    if (abs(imag(struct(i).alpha)) < sqrt(eps))
        struct(i).alpha = real(struct(i).alpha);
    end

    if ~b_dec
        str=[ str ' ' num2str(struct(i).alpha) '|' struct(i).bin '> +'];
    else
        str=[ str ' ' num2str(struct(i).alpha) '|' num2str(bin2dec(struct(i).bin)) '> +'];
    end
end

str=str(1:end-1);
