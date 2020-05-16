function [xOpt,tOpt,cOpt,startTimes,finishTimes,resSlack,actualResPath,travMat] = optDetermAN_multi_res_LP(w,iMat,K,rCosts,xLB,xUB,R,tPenalty,dueDt,cutsetCons,startTime)         
    
%this function solves the deterministic activity network (DAN) model for the specified project


%input vars:
%w          --> deterministic work content vector (certainty equivelents if finding Jensen LB)
%iMat       --> incidence matrix of the graph
%K          --> vector of model parameters indicating how duration is related to resource allocation d = w/(x^K)
%rCosts     --> vector of marginal resource costs for the activities 
%xLB        --> lower bounds on resource allocations
%xUB        --> upper bounds on resource allocations 
%R          --> the total units of resource available (used for aggregate constraint(s))
%tPenalty   --> cost of tardiness per time unit
%dueDt      --> due date of the project
%cutsetCons --> cell array of vectors where each vector is a udc to be constrained for total resource useage
%startTime  --> the start time of the project (t1 = startTime) 


%output vars:
%xOpt          --> vector of optimal resource allocations (x*) for the DAN model
%tOpt          --> vector of optimal node realization times (t*) for the DAN model
%cOpt          --> optimal cost secured for the DAN model
%startTimes    --> vector of activity start times under optimal policy (x*,t*)
%finishTimes   --> vector of activity completion times under optimal policy (x*,t*)
%seqA          --> vector giving the order of completion for the activities   
%seqN          --> vector giving the order of realization for the nodes
%resSlack      --> row vector giving the amount of surplus (unallocated) resource for each udc in cutsetCons
%actualResPath --> cell array of vectors w/ each vector a udc that the actual solution passed through
%travMat       --> matrix of 3 cols where first col is act num, second col is
%                  number of node activated (if any) by completion of act in col 1, and third is
%                  completion time - travMat gives traversal path in
%                  ascending order by activity completion time


%dependencies: 
%calls function findTraversal
%calls function findSurplusResource
%calls function findCutsetCons


%determine the size of the incidence matrix & convert to 2 col format
if size(iMat,2) > 2    %iMat is in long form so convert to short
    iM = convertImat(iMat);
else   %iMat is already in short form
    iM = iMat;
    iMat = convertImat(iMat);
end


%intialization
tic
numArcs = size(iM,1);
numNodes = max(iM(:,2));
numResources = size(K,2);
tOpt = zeros(1,numNodes);



%build a cell array indexed on nodes where cell i
%contains activities terminating in node i
termA = {[]};
inDegree = [0];
col2 = iM(:,2)';
a = [];
for n = 2:numNodes
    a = setdiff(ismember(col2,n) .* [1:numArcs] , 0);
    termA{n} = a;
    inDegree(n) = size(a,2);
end

%build a cell array indexed on nodes where cell i
%contains activities emminating from node i 
startA = {[]};
outDegree = [];
col1 = iM(:,1)';
a = [];
for n = 1:numNodes-1
    a = setdiff(ismember(col1,n) .* [1:numArcs] , 0);
    startA{n} = a;
    outDegree(n) = size(a,2);
end





%set the step size for discretization of the activity res allocations
lamdaStep = 0.5;


%construct the parameters for the LP (min cz st: Az=b)

%large number used for so that only 1 resource type is assigned to each
%activity (1 out of numResources binary integer constraints)
M = max(max(w))*1E6;

p = 0.4;
numVars = 0;
startPositions = zeros(1,numArcs);
endPositions = zeros(1,numArcs);
Aleft = [];
Aright = [];
A = [];
Aeqleft = [];
Aeqright = [];
AeqBottom = [];
Aeq = [];
c = [];
cLeft = [];
cRight = [];

for r = 1:numResources
    p = K(r);
    tempA = [];
    tempAeq = [];
    tempc = [];
    numVars = 0;
    for i = 1:numArcs
        xGrid = [xLB(r,i):lamdaStep:xUB(r,i)];
        numGridPoints(r) = size(xGrid,2);
        xGridInv = xGrid .^ (-1*p);
    %     xGridInv = xGrid .^ -1;
        temp = rCosts(r,i)*w(r,i)*xGrid.^(1-p);
        tempInv = w(r,i)*xGridInv;
        startIndex = numVars+1;
        startPositions(i) = startIndex;
        endIndex = numVars + size(temp,2);
        endPositions(i) = endIndex;
        tempc(startIndex:endIndex) = temp;
        tempA(i,startIndex:endIndex) = tempInv;
        tempAeq(i,startIndex:endIndex) = 1;   %cons so lamdas sum to 1
        numVars = numVars + size(temp,2);   
    end
    A = blkdiag(A,tempA);
    Aeq = blkdiag(Aeq,tempAeq);
    AeqBottom = [AeqBottom eye(numArcs)];
    c = [c tempc];
    Aleft = [Aleft;-M*eye(numArcs)];
    Aright =[Aright;iMat];
    
end

%append binary variables variables to A
A = [-M*eye(numArcs*numResources) A];
c = [zeros(1,numArcs*numResources) c];
Aeq = [[zeros(numArcs*numResources) Aeq];[AeqBottom zeros(numArcs,size(Aeq,2))]];

%append node realizations and tardiness vars
A = [A Aright zeros(size(A,1),1)];
c = [c zeros(1,numNodes) tPenalty];
Aeq = [Aeq zeros(size(Aeq,1),numNodes+1)];
% Aeq = [Aeq zeros(size(Aright)) zeros(size(Aeq,1),1)];

%add row for tardiness con
tCon = [zeros(1,size(A,2)-2) 1 -1];
A = [A;tCon];

%add cons for durations
% Ad = [eye(numArcs) zeros(numArcs,numResources*numGridPoints*numArcs) iMat zeros(numArcs,1)];
% A = [Ad;A];

%adjust variable counter for node realizations and tardinees
numVars = size(A,2);

%construct rhs vec for ineq cons
b = [zeros(1,size(A,1)-1) dueDt]';

%construct rhs vec for eq cons
beq = [ones(1,numArcs*numResources) (numResources-1)*ones(1,numArcs)]';
% beq = ones(1,size(Aeq,1))';

%add an eq con for node 1 (usually t1=0)
Aeq(size(Aeq,1)+1,numVars-numNodes) = 1;
beq(size(Aeq,1)) = startTime;

%establish the bounds on the vars
lb=zeros(1,numVars);
ub=[ones(1,numArcs*numResources) ones(1,numVars-(numArcs*numResources)-numNodes-1) startTime inf*ones(1,numNodes)];

%define vars that must be int valued
intcon = [1:numArcs*numResources];

%solve the linear program
[z,fval,exitflag,output,lambda] = linprog(c,A,b,Aeq,beq,lb,ub);
[zInt,fvalInt,exitflagInt,outputInt] = intlinprog(c,intcon,A,b,Aeq,beq,lb,ub);

%convert optimal solution from lamda form to resource allocations
startPositions = startPositions + numArcs*numResources;
endPositions = endPositions + numArcs*numResources;

endIndex = numArcs*numResources;

for r = 1:numResources
    for i = 1:numArcs
        startIndex = endIndex + 1;
        xGrid = [xLB(r,i):lamdaStep:xUB(r,i)];
%         startIndex = startPositions(i);
%         endIndex = endPositions(i);
        endIndex = startIndex + size(xGrid,2) - 1;
        xOpt(r,i) = sum(zInt(startIndex:endIndex)' .* xGrid);
    end
    if (r<numResources)
        startIndex = endIndex + 1;
%         startPositions = startPositions + numArcs*numGridPoints(r);
%         endPositions = endPositions + numArcs*numGridPoints(r);
    end
end

xOpt = vpa(xOpt);


%compute the optimal node realization times using Bellman Ford
tOpt = zeros(1,numNodes);
% startIndex = endPositions(numArcs) + 1; 
startIndex = endIndex + 1;
endIndex = numVars - 1;
tOpt(1:numNodes) = vpa(zInt(startIndex:endIndex));

cOpt = fval;

[tOptLP,cOptLP,criticalPathLP,piOptLP,z2,fval2,exitflag2,output2,lambda2] = linProgProjectEval(xOpt,w,iM,p,startA,termA,tPenalty,dueDt,startTime);
% [z,fval,exitflag,output,lambda] = linProgProjectEval(xOpt,w,iM,p,startA,termA,tPenalty,dueDt,startTime);

[tOptBF,cOptBF,criticalPathBF,piOptBF] = bellmanFordProjectEval(xOpt,w,iM,p,startA,termA,tPenalty,dueDt,startTime);

% dlmwrite('testA', A);
% dlmwrite('testAeq', Aeq);
% dlmwrite('testc', c);
% dlmwrite('testbeq', beq);
% dlmwrite('testb', b);


%need to add code to ensure that times that should be simultaneous are
%actually recognized as such            



% %check tOpt using Bellman Ford (structure of GP allows for alternate
% %optimal solns for tOpt for nodes off the critical path)
% [tOptBF,cOptBF,criticalPath,piOpt] = bellmanFordProjectEval(xOpt,w,iM,startA,termA,tPenalty,dueDt,startTime);
% 
% %    [tOptBF,cOptBF] = bellmanFordProjectEval(xOpt,w,iMat,rCosts,tPenalty,dueDt,startTime);
% timeDiff = tOpt - tOptBF;
% problemN = setdiff((abs(timeDiff) >= 1E-3*ones(size(timeDiff))) .* [1:numNodes] , 0);
% 
% %for any nodes that have different realization times than Bellman Ford
% %change to BF result - no delays in the start of seq feasible acts
% tOpt(problemN) = tOptBF(problemN);


%output activity start and completion times
d = w./(xOpt.^p);
s = iM(:,1)';
v = iM(:,2)';
startTimes = tOpt(s);
finishTimes = startTimes + d;

%round times off to approx 2 decimal places
startTimes = round(startTimes*50)/50;
finishTimes = round(finishTimes*50)/50;
tOpt = round(tOpt*50)/50;

%determine order of activity completion and node realization
seqA = sortrows([[1:numArcs]' finishTimes'],2);
seqA = seqA(:,1)';
seqN = sortrows([[1:numNodes]' tOpt'],2);
seqN = seqN(:,1)';

%output udc traversal path matrix
[travMat] = findTraversal(iM,tOpt,finishTimes,seqA,startTime);

%determine final udc traversal resource path 
[actualResPath] = findCutsetCons(iM,travMat);

%calculate the surplus resource for each udc on the path
[resSlack] = findSurplusResource(xOpt,actualResPath,R); 

%record computer run time to solve the problem
toc