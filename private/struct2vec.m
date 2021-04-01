function vec = struct2vec(struct)

n=length(struct(1).bin);
vec=zeros(2^n, 1);

for i=1:length(struct)
    dec = bin2dec(struct(i).bin);
    vec(dec+1)=struct(i).alpha;
end

