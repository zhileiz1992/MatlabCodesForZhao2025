parfor jj=200:209
%     parfor jj=1:48
      %   parfor jj=1:24
%       fs = flatnessAll(1).fs;
      % identify what files the onset and offset belong to
      % add padding, caution about boundary condition
      frameStart = max([1, onsets(warbleStart(jj))-padFramesBefore]);
      frameEnd = min([size(flatnessFlat,1), offsets(warbleEnd(jj))+padFramesAfter]);
      for fii=1:length(fileFrameIdx)
        if fii==1
          prev = 1;
        else
          prev = fileFrameIdx(fii-1);
        end
        if (frameStart>=prev) && (frameStart<=fileFrameIdx(fii))
          fileIdxStart = fii;
        end
        if (frameEnd>=prev) && (frameEnd<=fileFrameIdx(fii))
          fileIdxEnd = fii;
          break;
        end
      end
      fprintf('%d: %d, %d, %d\n', jj, fileIdxStart, fii, fileIdxEnd);
%       fileInRange = [fileStructAll(fileIdxStart:fileIdxEnd)];
end