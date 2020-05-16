function  [t,projectCost,criticalPath,piOpt] = bellmanFordProjectEval(xHat,w,Imat,K,startA,termA,tPenalty,T,t1)


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

        %find the longest path through the graph (makespan)
        %use shortest path alg with neg arcleng
        Et = zeros(1,numNodes);
        Et(1) = t1;
        for j=2:numNodes
            mx = -inf;
            a = termA{j};
            for i=1:size(a,2)
                k = Imat(a(i),1);
                test = Et(k) + d(a(i));
                if test > mx    
                    mx = test;
                end
            end
            Et(j) = mx;
        end

        %compute latest node realization times
        Lt = zeros(1,numNodes);
        Lt(numNodes) = Et(numNodes);
        for k = numNodes-1:-1:1
            minVal = inf;
            a = startA{k};
            for i = 1:size(a,2)
                j = Imat(a(i),2);
                test = Lt(j) - d(a(i));
                if test < minVal     
                    minVal = test;
                end
            end
            Lt(k) = minVal;
        end

        %return the early event times
        t(r,:) = Et;

        %determine if the project finished early
        tn = Et(numNodes);
        if tn - T <= zeroTol    %project is early
            projectCost(r) = 0;
        else     %project is tardy
            %determine the tardiness cost
            projectCost(r) = tPenalty * max([tn-T,0]);
            %return;  %exit function all dual vars = 0 
        end

        %determine the total floats for each activity
        TF = zeros(1,numArcs);
        for i = 1:numArcs
            a = [i Imat(i,:)];
            TF(i) = Lt(Imat(i,2)) - Et(Imat(i,1)) - d(i);
            if abs(TF(i)) <= zeroTol  %act is critical
                criticalPath = [criticalPath;a];
                piOpt(r,i) = -tPenalty;
            end
        end

    end
      
	
end