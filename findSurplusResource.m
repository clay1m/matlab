function [resSlack] = findSurplusResource(xOpt,cutsetCons,R)         
    
    %this function determines the total resource allocation for each udc and returns
    %the amount of surplus resource for each stage
    
    
    %input vars:
    %xOpt       --> vector of optimal resource allocations (x*) for the DAN model
    %cutsetCons --> cell array of vectors where each vector is a udc to be constrained for total resource useage
    %R          --> the total units of resource available (used for aggregate constraint(s))
    
    
    %output vars:
    %resSlack   --> row vector giving the amount of surplus (unallocated) resource for each udc in cutsetCons
    
    
    
    %determine the number of udcs on the resource useage path
    numU = max(size(cutsetCons));
    resSlack = zeros(1,numU);
    
    %calculate the surplus resource for each udc
    for u = 1:numU
        
        udc = cutsetCons{u};
        
        %calculate resource useage for the current udc
        totR = sum(xOpt(udc));
        
        %calculate the surplus resource
        resSlack(u) = R - totR;
        
    end
