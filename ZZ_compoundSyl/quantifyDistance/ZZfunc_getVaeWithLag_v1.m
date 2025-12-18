function [res, i_all] = ZZfunc_getVaeWithLag_v1(vae, win_loc, lags)
% get the vae latents with different lags to a list of spikes (win_loc)

res = nan(length(win_loc), length(lags), size(vae,2));

% construct a matrix for index
i_all = repmat(win_loc, 1, length(lags));
i_all = i_all + lags;

% if out of range, assign it to the first or last window
i_all(i_all<1) = 1;
i_all(i_all>size(vae,1)) = size(vae,1);

% grab values
res = vae(i_all,:);
res = reshape(res, [size(i_all,1) size(i_all,2) size(vae,2)]);


end

