noisy_pend = 'noisy_pendulum.csv';
true_pend =  'true_pendulum.csv';

f = noisy_pend;
Z = csvread(['data/' noisy_pend]);
Zperfect = csvread(['data/' true_pend]);
T = size(Z,1); % how many observations

name = sprintf('Estimate_noisy_%d.mp4',N);


writerObj = VideoWriter(name,'MPEG-4');
writerObj.FrameRate = 5;
% FrameRate
open(writerObj);


set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');


N = 500;
p = zeros(N,x_dims);

transitionprob = @(x0,noise) dynamic_model_pendulum(x0,noise,dt);
evalprob = @(z,x,R) evaluation_probability(z,x,R);


% generate initial particles from first observation
x0 = [Z(1,:), 0 , 0];
x0_noise = [0.2,0.2,0.7,0.7];
for i = 1:N
   p(i,:) = dynamic_model_pendulum(x0,x0_noise,1); 
end


%
% update the particles over time:
h = figure(1)
estimates = zeros(T,x_dims);
for t = 1:T
    
    clf
%     plot(p(:,1),p(:,2),'ro','MarkerSize',3); hold on
    plot(Z(t,1),Z(t,2),'b.','MarkerSize',10); hold on
    axis([-2 2 -4 1])
    
    p2 = zeros(N,x_dims);
    
    % motion update
    for n=1:N
       p2(n,:) = transitionprob(p(n,:),x0_noise);
    end
    p = p2;
    
    % plot the cloud
    plot(p(:,1),p(:,2),'b.','MarkerSize',2)
    
    % measurement update
    w = zeros(N,1);
    z_t = Z(t,:);
    for n = 1:N
       w(n) = evalprob(z_t,p(n,:),x0_noise(1:2));
    end
    % normalise:
    a = w/sum(w);
    
    % find max weighted particle and plot
    [val,i] = max(a);
%     plot(p(i,1),p(i,2),'rx','MarkerSize',10);
    
    % resampling (dumb algorithm)
    r1 = rand();
    
    
    % resampling (circular algorithm)
    p3 = zeros(N,x_dims);
    maxw = max(w(:));
    index = round(rand() * N); 
    % avoid getting 0 index 
    if (index == 0)
        % resample from 1 to N-1
        index = 1+round(rand() * N-1);
    end
    beta = 0;
    for n = 1:N
       beta = beta + (rand() * 2 * maxw);
       while beta > w(index)
          beta = beta - w(index);
          if (index == N-1)
              % then jump to start
              index = 1;
          else
              index = mod(index+1,N); %(never want to have 0);
          end
          
       end
       p3(n,:) = p(index,:);
        
    end
    
    p = p3;
    
    % plot the cloud
    plot(p(:,1),p(:,2),'ro','MarkerSize',2)
    % plot the best estimate (mean)
    % maybe should use weighted mean?

    best_estimate = mean(p,1);
    estimates(t,:) = best_estimate;
    plot(best_estimate(1),best_estimate(2),'g+','MarkerSize',15)
    hold off
    legend('measurement','old particles','resampled particles','estimate');
%     pause(0.1)
    frame = getframe;
    writeVideo(writerObj,frame);

    
end


close(writerObj);