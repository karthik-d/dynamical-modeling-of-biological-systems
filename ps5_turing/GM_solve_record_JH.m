%%    GM_solve_record_JH
%     
%in this script, we write a function that solves the Gierer_Meinhardt 
%system for parameters D, omega and sigma. We use MATLAB's PDEPE to do so.
%The function also plots the solution 1) as a surface & 2) as an animation.
%The function then saves the animation also.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol] = GM_solve_record_JH(P,tmax,delt,L,filename)
    %we solve the PDE for parameters
    %P(1) = D in PDE;
    %P(2) = omega in PDE;
    %P(3) = sigma in PDE;
    %tmax = amount of time that we run the system;
    %L = length (spatial) of system
    %filename = name of file to save the animation to, e.g. 'patterns.avi'
    
    GMfunctions = GMfuns_JH; %here we call the necessary functions to solve the PDE:
                     %the initial conditions function, the boundary conditions
                     %function and the PDE function itself.
    rng(404); %set the seed for random number generation
    t = linspace(0,tmax,tmax/delt); %time vector (n time points)
    x = linspace(0,L,200); %the mesh on which we solve (200 space points)
    %%PDEPE solver
    sol = pdepe(0,@GMfunctions.fun3,@GMfunctions.fun1,@GMfunctions.fun2,x,t,[],P);
    u1 = sol(:,:,1); %activator solution
    u2 = sol(:,:,2); %inhibitor solution



    %% Plotting in 3 DIMENSIONS

    figure(1)
    surf(x,t,u1,'edgecolor','none');
    xlabel('Position','fontsize',20,'fontweight','b','fontname','arial')
    ylabel('Time','fontsize',20,'fontweight','b','fontname','arial')
    zlabel('[Activator]','fontsize',20,'fontweight','b','fontname','arial')
    axis([0 L 0 tmax 0 max(max(u1))])
    set(gcf(), 'Renderer', 'painters')
    set(gca,'FontSize',18,'fontweight','b','fontname','arial')

    figure(2)
    surf(x,t,u2,'edgecolor','none');
    xlabel('Position','fontsize',20,'fontweight','b','fontname','arial')
    ylabel('Time','fontsize',20,'fontweight','b','fontname','arial')
    zlabel('[Inhibitor]','fontsize',20,'fontweight','b','fontname','arial')
    axis([0 L 0 tmax 0 max(max(u2))])
    set(gcf(), 'Renderer', 'painters')
    set(gca,'FontSize',18,'fontweight','b','fontname','arial')

    %%
    writerObj = VideoWriter(filename);
    writerObj.FrameRate = 5;
    open(writerObj);
    figure(3)
    for n = 2:2:length(t)
        set(gca, 'FontSize', 18, 'LineWidth', 1); %<- Set properties
        plot( x , sol(n,:,1), 'LineWidth',3);
        hold on
        plot( x , sol(n,:,2), 'r', 'LineWidth',3);
        hold off
        legend('activator', 'inhibitor', 'Location', 'SouthEast');
        title(strcat('Gierer-Meinhardt patterns t =' , sprintf(' %d ', ceil(t(n)))));
        axis([0 L 0 max(max(max(sol(:,:,:))))+0.1])
        M(n) = getframe(figure(3));
        writeVideo(writerObj,M(n));
    end

    % play as smooth movie 1 time at 5 frames per second
%     numtimes=1;
%     fps=5;
%     movie(M,numtimes,fps)
end