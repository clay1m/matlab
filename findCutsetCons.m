function [cutsetCons] = findCutsetCons(iMat,travMat)         
    
    %this function determines the builds the cutset constraint cell array corresponding to a particular
    %udc traversal path for a given project realization
    
    
    %input vars:
    %iMat       --> incidence matrix of the graph
    %travMat    --> matrix of 3 cols where first col is act num, second col is
    %               number of node activated (if any) by completion of act in col 1, and third is
    %               completion time - travMat gives traversal path in
    %               ascending order by activity completion time
    
    %output vars:
    %cutsetCons --> cell array of vectors where each vector is a udc to be constrained for total resource useage
    
    
    %initialize vars
    numArcs = size(iMat,1);
    numNodes = max(iMat(:,2));
    seqA = travMat(:,1)';
    nodeCol = travMat(:,2)';
    timeCol = travMat(:,3)';
    s = iMat(:,1)';
    v = iMat(:,2)';
    cutsetCons = {};
    udcCount = 1;
    U = setdiff(ismember(s,1) .* [1:numArcs] , 0);
    cutsetCons{udcCount} = U;
    udcCount = udcCount + 1;
    
    
    %build a list of of distinct event times
    
    eventTimes = union(timeCol,timeCol);  %remove duplicate times from event list
    
    for e = 1:size(eventTimes,2)
        
        %mark the time of the current event
        t = eventTimes(e);
        
        %determine the rows of travMat that are invlved in the current event
        I = setdiff(ismember(timeCol,t) .* [1:numArcs] , 0);
        
        %determine act(s) completed
        cA = seqA(I);
        
        %remove the completed act(s) from the active udc
        U = setdiff(U,cA);
        
        %determine if any nodes were realized
        rN = setdiff(nodeCol(I) , 0);
        
        %if node(s) are realized add new sequence feasible acts to the active udc
        Asf = [];    %new sequence feasible acts
        for j = 1:size(rN,2)
            n = rN(j);
            Asf = union(Asf , setdiff(ismember(s,n) .* [1:numArcs] , 0));
        end
        U = union(U,Asf);
        
        %record the active udc
        if size(U,2) > 0       %do not record empty sets
            cutsetCons{udcCount} = U;
            udcCount = udcCount + 1;
        end
        
    end
    celldisp(cutsetCons)
    
    