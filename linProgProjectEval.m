function  [t,projectCost,criticalPath,piOpt,z,fval,exitflag,output,lambda] = linProgProjectEval(xHat,w,Imat,K,startA,termA,tPenalty,T,t1)

	%this function determines the node realization times and the tardiness 
    %cost of the project for a given workcontent vector w and a resource 
    %allocation x using the Bellman Ford shortest path algorithm with 
    %negative arc lengths. Using crtical path info this function generates
    %the simplex multipliers for use in Bender's decomp w/o solving the LP
	
	
	%input vars:
	%xHat     --> vector of resource allocations for each activity
	%w        --> deterministic work content vector (or matrix)
	%Imat     --> incidence matrix of the graph (in 2 col short form)
    %K        --> model parameter indicating how duration is related to resource allocation d = w/(x^K)
    %startA   --> cell array - cell i contains set of acts starting in node i
    %termA    --> cell array - cell i contains set of acts ending in node i
	%tPenalty --> cost of tardiness per time unit
	%T        --> due date of the project
	%t1       --> the start time of the project
	
	
	%output vars:
	%t             --> vector of node realization times for the DAN model
	%projectCost   --> project cost for the DAN model
    %criticalPath  --> list of activities on the critical path (or subgraph)
	%piOpt         --> dual variables from LP formulation
    %z             --> LP soln for node realizations times
    %fval          --> LP optimal function value (project completion time)
    %exitflag      --> LP algorithm exit condition
    %output        --> LP algorithm output
    %lambda        --> LP algorithm reduced costs
	
	%set the tolerance for zero values
    zeroTol = min([mean(w)/100 1E-2]);
    
    %check if format of imat is correct
    if size(Imat,2)>2
        Imat = convertImat(Imat);
    end
    
	%intialization
	numArcs = size(Imat,1);     
	numNodes = max(Imat(:,2));
    numRuns = size(w,1);
    criticalPath = [];
    t = zeros(numRuns,numNodes);
    projectCost = zeros(1,numRuns);
    Kt = zeros(1,numRuns);
    piOpt = zeros(numRuns,numArcs);
    
    
     
    %loop over thr number of w vectors input by user
    
    for r = 1:numRuns
	    
        %show iteration count
        r
        
        %calculate the durations of the activities
        d = w(r,:)./(xHat.^K);
        
        %construct the parameters for the LP (min cz st: Az=b)
        A = zeros(numArcs,numNodes);
        b = zeros(numArcs,1);
        c = [zeros(1,numNodes-1) 1];
        c(1,1) = -1;
       
        for i = 1:numArcs
            arc = Imat(i,:);
            startNode = arc(1);
            endNode = arc(2);
            A(i,startNode) = 1;
            A(i,endNode) = -1;
            b(i,1) = -d(i);
        end

        %adjust variable counter for node realizations and tardinees
        numVars = size(A,2);

        
        %add an eq con for node 1 (usually t1=0)
%         Aeq(numArcs+1,numVars-numNodes-1) = 1;
%         beq(numArcs+1) = startTime;

        %add cols of zeros to make Aeq the correct dimensions
%         Aeq = [Aeq zeros(numArcs+1,numNodes+1)];

        %establish the bounds on the vars
%         lb=-inf*ones(1,numVars);
        lb = zeros(1,numVars);
        ub= inf*ones(1,numVars);

        %solve the linear program
        [z,fval,exitflag,output,lambda] = linprog(c,A,b,[],[],lb,[]);
    end
	
end