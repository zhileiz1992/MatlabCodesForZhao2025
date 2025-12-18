function [is_pass] = ZZfunc_calcRange_v1(linear_i, size_d, min_dur)

  % convert to xy index
  [xx, yy] = ind2sub(size_d, linear_i);
  x_range = max(xx) - min(xx); 
  y_range = max(yy) - min(yy); 
  is_pass = (x_range>=min_dur) & (y_range>=min_dur);

end

