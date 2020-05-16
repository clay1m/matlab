function [ImatOut] = convertImat(ImatIn)

%This function converts between long and short forms of the incidence
%matrix
%
%INPUT
%ImatIn    ---> Incidence matrix for the graph
%
%OUTPUT
%ImatOut   ---> Incidence matrix for the graph converted to alaternate form

%determine the number of arcs
numArcs = size(ImatIn,1);     

%determine which form the incidence matrix currently is
test = size(ImatIn,2);
if test == 2   %short form so convert to long
	numNodes = max(ImatIn(:,2));
    ImatOut = zeros(numArcs,numNodes);
    for i = 1:numArcs
        a = ImatIn(i,:);
        ImatOut(i,a(1)) = 1;
        ImatOut(i,a(2)) = -1;
    end
else   %long format so convert to short
    
    ImatOut = zeros(numArcs,2);
    numNodes = size(ImatIn,2);
    test = ones(1,numArcs)' * [1:numNodes];
    temp = ImatIn ~= zeros(size(ImatIn));
    temp = temp .* test;
    
    for i = 1:numArcs
        ImatOut(i,:) = setdiff(temp(i,:),0);
    end
    
end


    
    
    


end

