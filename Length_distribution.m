clear all;
clc

er=1e-4;
nm=370;

fr="5";
cp="5";

c=load(['L_data/L_fr_' char(fr) '_cp_' char(cp) '.txt']);
c(c==0) = [];
l=c/nm;

n=max(l);
maxL = max(l);
binw = 0.5;


%%%%%%%%%% P(n) %%%%%%%%%

pp=[];

nmax = n

%----------------------------------

binpos=[0:binw:nmax];
bincen=[binw/2:binw:nmax-binw/2];

[counts] = histcounts(l,binpos,'Normalization', 'pdf');	

pp = counts';

fid21=fopen(['data/dist' '_fr_' char(fr) '_cp_' char(cp)  '.txt'],'w');

for i=1:min((maxL/binw)+20,nmax/binw)
fprintf(fid21, '%f %f\n', bincen(i), pp(i));
end

fclose(fid21)

psum = sum(pp*binw)

figure(1)
plot(bincen, pp)












