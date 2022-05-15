clear all;
clc

er=1e-4;
nm=370;			% subunit to micro-meter conversion factor
usermax=300;		% time to be extracted
dt=5;			% time range of local average- avg over t to t+dt

fr="2";
cp="0";

c=load(['L_data/P_fr_' char(fr) '_cp_' char(cp) '.txt']);

%% time, number of filament, mean length in # , capping bound filaments #, profilin, profilin-actin, formin-bound filaments #, reaction count, monomer count

l=c(:,3)/nm;		% extract length
t=c(:,1);		% extract time

tmax=min(usermax,max(t));
tv=[0:dt:tmax];

tt=[];
ll=[];
lsig=[];

for i=1:numel(tv)-1
tt(i)=(tv(i)+tv(i+1))/2;
ls=[];
lsum=0;
lcnt=0;
for j=1:numel(l)
if (t(j)>=tv(i)&&t(j)<=tv(i+1))
lsum=lsum+l(j);
lcnt=lcnt+1;
ls=[ls, l(j)];
end
end
ll(i)=lsum/lcnt;
lsig(i)=std(ls)/sqrt(numel(ls));
end

% set 0 at time=0
tt=[0, tt];
ll=[0, ll];
lsig =[0, lsig];

fid21=fopen(['data/mean' '_fr_' char(fr) '_cp_' char(cp)  '.txt'],'w');

for i=1:numel(tt)
fprintf(fid21, '%f %f %f\n', tt(i), ll(i), lsig(i));
end

fclose(fid21)

errorbar(tt,ll,lsig)


















