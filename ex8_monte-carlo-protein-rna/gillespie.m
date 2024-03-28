function [outX,outT] = ...
    gillespie(X0,calculatePropensities,stochiometry,p)
%gillespie standard gillespie solver
%   p                          structure that contains the fixed parameters
%                              following entries must be defined
%                               p.startTime = start time
%                               p.endTime = end time
%                               p.nrep = number of replicate 
%                               p.dtForOutput = output time interval
% 
%   X0                         is a vector of initial condition
%
%   calculatePropensities(X,p) is a function that calculate propensities
%                              note if using additional constants, they
%                              need to be defined in p.
%
%   stochiometry               stochiometry matrix, 
%                               # rows = # species
%                               # columns = # reactions
%   
%
%
% check input and throw error if not ok
sz = size(X0);
if sz(1) < sz(2)
    error('initial condition X0 must be a coummn vector')
end
if sum(abs(X0-round(X0))) ~= 0
    error('entries in the initial condition X0 must be integers')
end
if sum(X0<0)>0
    error('entries in the initial condition X0 must be >=0 ')
end
X0=reshape(X0,[],1);%make a column vector
% create array to hold variables
X = zeros(length(X0),1);

% create arrays for output
outX = zeros(p.endTime/p.dtForOutput,p.nrep,length(X));
outT = zeros(p.endTime/p.dtForOutput,1); %aray to store time

   % generate our random number ahead of time for efficiency, save it
    % to a matrix called randNums.  We also need a counter to determine
    % which random number to use: Gillespie's algorithm uses 2 random numbers
    % at each step.
    numberOfRandomNums = 5000000;
    randNums = rand(numberOfRandomNums,2);
    rIndex = 1;


% run the simulation
for rep=1:p.nrep
    
    %initial condition
    time = p.startTime;
    X = X0;
    
    % bookeeping output
    nextOutputTime = time+p.dtForOutput; % set the time to output next
    outX(1,rep,:) = X;
    ouT(1) = time;
    lastResults = X;
    resultIndex = 2;
    
    %define the list of reaction propensities based on the initial
    %conditions and parameters.
    propensities = calculatePropensities(X,p);

    if sum(propensities < 0)>0
        error(['One of the propensities is negative. This can happen', ...
            'when you make use of conservation laws in the propensities ',...
            'and define the total number of one of the species using a ',...
            'real number rather than an integer. Make sure that you define ',...
            'total numbers of molecules as integers. You can use ',...
            'round(...) to make sure the number is an integer.'])
    end
    
 
    
    % begin the main simulation loop of the algorithm
    while time <= p.endTime
        
        % step 0, if we used up all random numbers, generate a new batch
        if ( rIndex > numberOfRandomNums )
            % fprintf('time = %4.2f Used all random numbers. Generating new batch\n',time);
            randNums = rand(numberOfRandomNums,2);
            rIndex = 1;
        end
        
        % step 1, generate the time to the next event (requires us
        % to calculate the total propensity of the system! The total
        % propensity is used to generate a random number exponentialy
        % distributed).
        rtot = sum(propensities);
        tau = -(1/rtot)*log(randNums(rIndex,1));
        
        % step 2, select the next reaction to fire- first we have to
        % draw a random number between 0 and 1 the total propensity sum
        %value = rand(1)*rtot; %
        value = randNums(rIndex,2)*rtot;
        
        % step2 / step 3, actually fire the reaction by updating the
        % counts based on 'value' that we generated.
        
        %find out which reaction fires
        reaction = sum(cumsum(propensities) < value)+1;
        
        %execute that reaction
        X = X + stochiometry(:,reaction);
        
        
        % step 4, update the propensities and the simulation time
        propensities = calculatePropensities(X,p);
        if sum(propensities < 0)>0
            error(['One of the propensities is negative. This can happen', ...
                'when you make use of conservation laws in the propensities ',...
                'and define the total number of one of the species using a ',...
                'real number rather than an integer. Make sure that you define ',...
                'total numbers of molecules as integers. You can use ',...
                'round(...) to make sure the number is an integer.'])
        end
        time = time + tau;
        rIndex = rIndex+1;
        
        % step 5, generate results from the simulation by saving it to
        % our results output matrix.  Only generate output when the simulation
        % time exceeds the nextOutputTime.
        if time>nextOutputTime
            
            % if we jumped over the nextOutputTime we need to update our array
            while time>nextOutputTime
                outX(resultIndex,rep,:) = lastResults;
                outT(resultIndex) = nextOutputTime;
                nextOutputTime = nextOutputTime + p.dtForOutput;
                resultIndex = resultIndex+1;
                if nextOutputTime>p.endTime
                    break;
                end
            end
            
        end
        
        % store results of last jump in case next we jump beyond storage point
        % and need to store intermediate time points
        lastResults = X;
        
    end
    
end

%fix the output
outT = outT(1:end-1);
outX = outX(1:end-1,:,:);

end

