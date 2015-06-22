function r_vec_out = brownian_motion(time, delta_x, dtstep, rmax, walkers)
%function [r_vec_out,z_vec_out] = brownian_motion(time,delta_t,rmax, walkers)

%This function spits out the radius of the position of a number of walkers
%who are dropped randomly in a circle of radius rmax and move around
%randomly in the circle over time. The output is a vector of size Nstep, 
% where Nstep = time/dt, where dt=dtheta^2 rmax^2. Thus the ratio 
% T/rmax^2 = Nstep*dtheta^2. This number governs if the randomwalk has
% converged or not. 
% typical values for dtheta=0.01 or 0.03

%time = 100;
%dtheta = 0.03;
%rmax = 1;
%dtstep = rmax^2*dtheta^2; %counts steps of output.
%walkers = 1;

nsteps = round(time/dtstep);

%msteps = round(dtstep/rmax^2/dtheta^2);
%This governs the inner loop of the walk

%r_vec_out = zeros(nsteps,walkers);
x_vec_out = zeros(nsteps,walkers);
y_vec_out = zeros(nsteps,walkers);


%initialize positions
r_vec_init = rmax*sqrt(rand(1,walkers));
theta_vec_init = 2*pi*rand(1,walkers);


%r_vec_out(1,:) = r_vec_init;
x_vec_out(1,:) = r_vec_init.*cos(theta_vec_init);
y_vec_out(1,:) = r_vec_init.*sin(theta_vec_init);

%this figure shows distribution of initial users.
%figure(3); plot(x_vec_out(1,:), y_vec_out(1,:),'.');

%outer loop
for ii = 1:(nsteps-1)
    %inner loop
    
    xtemp = x_vec_out(ii,:);
    ytemp = y_vec_out(ii,:);
%    for jj = 1:msteps
        
        step_xvec = delta_x*(-3+2*randi(2,[1,walkers]))/sqrt(2);
        xtemp = xtemp + step_xvec;
        step_yvec = delta_x*(-3+2*randi(2,[1,walkers]))/sqrt(2);
        ytemp = ytemp + step_yvec;
        
        ind_vec = find(abs(xtemp+1i*ytemp)>rmax);
        if ~isempty(ind_vec)
%             xtemp1 = xtemp(ind_vec) - 0.5*step_xvec(ind_vec);
%             ytemp1 = ytemp(ind_vec) - 0.5*step_yvec(ind_vec);
%             theta1 = atan2(ytemp1,xtemp1);
%             xtemp(ind_vec) = xtemp1 + 0.5*step_xvec(ind_vec) - delta_x*cos(theta1);
%             ytemp(ind_vec) = ytemp1 + 0.5*step_yvec(ind_vec) - delta_x*sin(theta1);

            xtemp(ind_vec) = xtemp(ind_vec) - step_xvec(ind_vec);
            ytemp(ind_vec) = ytemp(ind_vec) - step_yvec(ind_vec);
            if size(find(abs(xtemp+1i*ytemp)>rmax)>0)
                qqq = 'There is a problem';
                return
            end
        end
        
%    end
    x_vec_out(ii+1,:) = xtemp;
    y_vec_out(ii+1,:) = ytemp;
end

r_vec_out = abs(x_vec_out+1i*y_vec_out);

circ_vec = rmax*[cos(2*pi*(0:nsteps)'/nsteps) sin(2*pi*(0:nsteps)'/nsteps)];

figure(1); clf; hold on;
% for ii = 1:walkers
% plot(x_vec_out(:,ii),y_vec_out(:,ii),'g', circ_vec(:,1), circ_vec(:,2),'r' )
% end


plot(x_vec_out(:,1),y_vec_out(:,1),'k', circ_vec(:,1), circ_vec(:,2),'k' )
plot(x_vec_out(:,2),y_vec_out(:,2),'k', circ_vec(:,1), circ_vec(:,2),'k' )
plot(x_vec_out(:,3),y_vec_out(:,3),'y', circ_vec(:,1), circ_vec(:,2),'k' )
plot(x_vec_out(:,4),y_vec_out(:,4),'c', circ_vec(:,1), circ_vec(:,2),'k' )
plot(x_vec_out(:,5),y_vec_out(:,5),'r', circ_vec(:,1), circ_vec(:,2),'k' )
plot(x_vec_out(:,6),y_vec_out(:,6),'g', circ_vec(:,1), circ_vec(:,2),'k' )

figure(3); clf;
plot(x_vec_out(:,1),y_vec_out(:,1),'k', circ_vec(:,1), circ_vec(:,2),'k' )
% figure(2); plot(1:nsteps,r_vec_out(:,1))
