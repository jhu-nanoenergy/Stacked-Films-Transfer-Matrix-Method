function [outputStruct,thickness]=ParticleSwarmOptm(obj)
f=@(t) ColorDelta(t,obj);
ub=[300 300 600 150 200];
lb=[ 50 50 250 10 100];
options=optimoptions('particleswarm','SwarmSize',70,'Display','iter','FunctionTolerance',1e-5,'displayinterval',3,'UseParallel',1);
%Options of the particle swarm algorithm with 140 particles, displaying iterative results after each 3 iterations, using parallel computing version
tic
[thickness,f,flag,output]=particleswarm(f,5,lb,ub,options);
outputStruct=Objective6layersPbI2_Au(thickness);
% x is the optimal thickness
toc