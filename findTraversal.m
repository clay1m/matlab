function [travMat] = findTraversal(iMat,tOpt,finishTimes,seqA,startTime)         
    
    %this function determines the udc traversal path for a given project realization
    
    
    %input vars:
    %iMat       --> incidence matrix of the graph
    %tOpt       --> vector of optimal node realization times (t*) for the DAN model
    %finishTimes --> vector of activity completion times under optimal policy (x*,t*)
    %seqA       --> vector giving the order of completion for the activities 
    %startTime  --> the start time of the project (t1 = startTime) 
    
    %output vars:
    %travMat    --> matrix of 3 cols where first col is act num, second col is
    %               number of node activated (if any) by completion of act in col 1, and third is
    %               completion time - travMat gives traversal path in
    %               ascending order by activity completion time
 




    %initialize vars
    timeCol = finishTimes(seqA);
    nodeCol = [];
    numArcs = size(iMat,1);
    numNodes = max(iMat(:,2));
    v = iMat(:,2)';
    
    %loop through the activities to determine the traversal order
    for i = 1:numArcs
        
        a = seqA(i);
        t = timeCol(i);     %completion time of activity a
        rN = [];            %set of nodes realized at time t
        
        if t > startTime         %activity has positive duration (not a dummy)
            
            test = sum(ismember(tOpt,t));
            rN = setdiff(ismember(tOpt,t) .* [1:numNodes] , 0);
            
            if test == 0   %the act in question does not activate a node
                
                nodeCol(i) = 0;
                
            else           %activity does terminate at the realization time of one or more nodes
                
                if sum(ismember(rN,v(a))) > 0  %activity a triggers its terminal node
                    
                    nodeCol(i) = v(a);
               
                else       %activity completion does not activate a node
                    
                    nodeCol(i) = 0;
                    
                end
                
            end
        else                     %activity has zero duration (a dummy)
            nodeCol(i) = 0;
            seqA(i) = 0;
            
        end
        
    end
    
    travMat = [seqA' nodeCol' timeCol'];  
    travMat = sortrows(travMat,[3 1]);