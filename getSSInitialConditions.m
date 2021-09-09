function init=getSSInitialConditions(preInit,v,absThresholdForSS,timeStepForSS)
    
    v.drugTC=ones(1,timeStepForSS);
    [t,x]=ode15s(@(t,x) mmODE(t,x,[],v), [0 timeStepForSS], preInit,v);
    
    while(max(abs(x(2,:)-x(end,:)))>absThresholdForSS)
        fprintf('diff: %d, looping again to find SS\n',max(abs((x(2,:)-x(end,:)))));
        newInit=x(end,:);
        [t,x]=ode15s(@(t,x) mmODE(t,x,[],v), [0 timeStepForSS], newInit,v);
    end
    fprintf('diff:%d, ss found\n',max(abs((x(2,:)-x(end,:)))));
    init=x(end,:);
end