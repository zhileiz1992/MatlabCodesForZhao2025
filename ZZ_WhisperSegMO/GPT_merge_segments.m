function [t_merged, n_merged] = GPT_merge_segments(t, n, thre_pt)
    % Initialize
    t_merged = [];
    n_merged = {};
    
    % Start with the first segment
    curr_start = t(1,1);
    curr_end   = t(1,2);
    curr_name  = n{1};
    
    for i = 2:size(t,1)
        next_start = t(i,1);
        next_end   = t(i,2);
        next_name  = n{i};

        if next_start - curr_end < thre_pt
            % Merge segments
            curr_end = next_end;
            if ~strcmp(curr_name, next_name)
                curr_name = 'x';
            end
        else
            % Save the current segment
            t_merged = [t_merged; curr_start, curr_end];
            n_merged = [n_merged; {curr_name}];
            % Start a new segment
            curr_start = next_start;
            curr_end = next_end;
            curr_name = next_name;
        end
    end

    % Append the final segment
    t_merged = [t_merged; curr_start, curr_end];
    n_merged = [n_merged; {curr_name}];
end
