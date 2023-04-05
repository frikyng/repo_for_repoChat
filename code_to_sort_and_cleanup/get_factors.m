function LoadingsPM = get_factors(data)

    valid = ~all(isnan(data));
    data = fillmissing(data,'linear');
    
    [LoadingsPM,specVarPM,rotationPM,stats, scores] = ...
                    factoran(data(:, valid),10);
    
    step = floor(prctile(scores(:),99.9));            
    %figure();plot(scores - (1:step:size(scores, 2)*step))            
    
    
%     data_fa = scores*LoadingsPM'; 
%     data_fa = bsxfun(@times, data_fa, std(data));
%     data_fa = bsxfun(@plus, data_fa, mean(data)); 
%     
%     figure(123);
%     for trace = 1:size(data, 2)
%         figure(123);cla; hold on; plot(data(:, trace),'r', 'LineWidth',2); hold on;plot(data_fa(:, trace),'k');
%         drawnow;pause(0.1);
%     end
    
%     data_sub = data - data_fa(:,1);
%     
%     figure(123);
%     for trace = 1:size(data, 2)
%         figure(123);cla; hold on; plot(data(:, trace),'r', 'LineWidth',2); hold on;plot(data_sub(:, trace),'k');
%         drawnow;pause(0.1);
%     end
    
    
%     n_gp = size(data, 2);
% 
%     [W, H] = nnmf(data, n_gp,'replicates', 20);
%     %[W1, H1, T1] = factoran(data, 20);
%     [B,L,var,fac,E] = FA(data);
%     figure();plot(mean(W*H, 2) - mean(data, 2)); hold on; plot(mean(data, 2));

%     %% Remove temporal components that are to faint
%     keep = max(W) > max(W(:)/20);
%     W    = W(:, keep);
%     H    = H(keep, :);
%     n_gp = sum(keep);
%     
    %% Remove spatial component that are not widespread
    %keep = sum(H > 0.05, 2) > round(size(data, 2)/10);
    %W    = W(:, keep);
    %H    = H(keep, :);
    %n_gp = sum(keep);
   
    %hold on; plot(mean(W*H, 2), 'k');
    
%     figure();
%     axes = [];
%     for i = 1:n_gp
%         ax = subplot(n_gp,2,(i-1)*2+1);
%         plot(W(:,i));
%         subplot(n_gp,2,i*2);
%         plot(H(i,:));
%         axes = [axes, ax];
%     end
%     linkaxes(axes, 'x')

    
%     axes = [];
%     for i = 1:n_gp
%         ax = subplot(n_gp,2,(i-1)*2+1);
%         plot(W(:,i));
%         subplot(n_gp,2,i*2);
%         plot(H(i,:));
%         axes = [axes, ax];
%     end
%     linkaxes(axes, 'x')


end


