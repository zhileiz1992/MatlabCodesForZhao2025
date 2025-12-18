function [res, i_all] = ZZfunc_getVaeWithLag_v2(vae, win_loc, lags)
% get the vae latents with different lags to a list of spikes (win_loc)
% differ from v1: assign to nan if out of range

res = nan(length(win_loc), length(lags), size(vae,2));

% construct a matrix for index
i_all = repmat(win_loc, 1, length(lags));
i_all = i_all + lags;

% if out of range, assign it to the first or last window
for ii=1:size(i_all,1)
  for jj=1:size(i_all,2)
    if (i_all(ii,jj)>=1) && (i_all(ii,jj)<=size(vae,1))
      res(ii,jj,:) = vae(i_all(ii,jj),:);
    end
  end
end

end

