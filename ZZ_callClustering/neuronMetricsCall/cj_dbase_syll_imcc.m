function [cc] = cj_dbase_syll_imcc(dbase,syll,binsz)
    % ifr params
    gaussian_kernel = @(kt, kt_i, sigma) (1 / (sigma * sqrt(2 * pi))) * exp(-(kt - kt_i).^2 / (2 * sigma^2));
    sigma = 0.02; % 20 ms width gaussian

    fs = dbase.Fs;
    voc = dbase.voc_segs_zz;
    labs = dbase.whisper_seg_labels_zz;
    sortf = dbase.sortedfiles;
    sortedspks = dbase.sortedspks;

    count = 0;
    for i = 1:length(sortf)
        voci = find(voc.ifile==sortf(i));
        if isempty(voci)
            continue
        end
        theselabs = labs{voci};
        thesesylls = voc.whispersegs{voci};
        sylli = find(cellfun(@(x) isequal(x,syll),theselabs));
        % spks
        filespks = sortedspks{i};
        kt = 0:binsz/1000:20; % files are always 20s
        ifr = zeros(size(kt));
        for q = 1:length(filespks)
            ifr = ifr+gaussian_kernel(kt,filespks(q)/fs,sigma);
        end
        for k = 1:length(sylli)
            count = count +1;
            sylltime = thesesylls(sylli(k),:)/fs;
            durs(count) = sylltime(2)-sylltime(1);
            allifrs{count,:} = ifr(kt>sylltime(1)&kt<sylltime(2));
        end
    end

    % time warp ifrs to median duration
    meddur = median(durs);
    warped = cj_timewarp_vectors(allifrs,round(meddur/(binsz/1000)));

    % get mean of corr matrix
    R = corrcoef(warped');
    flat = reshape(R,[1,size(R,2)^2]);
    cc = mean(flat(~isnan(flat)));

end