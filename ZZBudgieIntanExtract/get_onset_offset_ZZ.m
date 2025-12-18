function [onsets,offsets]=get_onset_offset_ZZ(signal,loc,thres,gapsize)
% modified from Han's function
% allow a small gap that below the threshold when finding onset/offset
    onsets=[]; offsets=[];
    loc(loc==1)=2;
    loc(loc==length(signal))=length(signal)-1;
    onset=loc(1)-1;
    offset=loc(1)+1;
    while onset>1
        % check if no points pass the threshold within the window
        start_idx = max(1, onset-gapsize);
        check_pass = sum(signal(start_idx:onset)>thres);
        if check_pass>0
            onset=onset-1;
        else
            break;
        end
    end
    onsets =[onsets onset];
    while offset<length(signal)
        % check if no points pass the threshold within the window
        end_idx = min(length(signal), gapsize+offset);
        check_pass = sum(signal(offset:end_idx)>thres);
        if check_pass>0
            offset=offset+1;
        else
            break;
        end
    end
    offsets = [offsets offset];
    for i=2:length(loc)
        if  loc(i)>max(offsets)  
            onset=loc(i)-1;
            offset=loc(i)+1;
            while onset>1
                % check if no points pass the threshold within the window
                start_idx = max(1, onset-gapsize);
                check_pass = sum(signal(start_idx:onset)>thres);
                if check_pass>0
                    onset=onset-1;
                else
                    break;
                end
            end
            onsets =[onsets onset];
            while offset<length(signal)
                % check if no points pass the threshold within the window
                end_idx = min(length(signal), gapsize+offset);
                check_pass = sum(signal(offset:end_idx)>thres);
                if check_pass>0
                    offset=offset+1;
                else
                    break;
                end
            end
            offsets = [offsets offset];
        end
    end
end

