%% plot the call spectrogram
% use Brian's function
function [fig_count]=plotSpectrogramListBrian9Unscale(sound_struct, save_folder, expID, suffix, fig_start, maxDur, fs, plot_row, plot_col, width, height)
    % pad zero the short segments
    numPoint = floor(maxDur*fs);
    sound_structPad = sound_struct;
    for ii=1:length(sound_structPad)
        s = sound_structPad(ii).signalNorm;
        sPad = zeros(numPoint, 1);
        idxStart = max(1, floor(numPoint/2 - length(s)/2));
        idxEnd = min(numPoint, idxStart+length(s)-1);
        sPad(idxStart:idxEnd) = s(1:(idxEnd-idxStart+1));
        sound_structPad(ii).signalNorm = sPad;
        sound_structPad(ii).onsetTimeOfFile = 0;
        sound_structPad(ii).offsetTimeOfFile = maxDur;
    end
    sound_struct = sound_structPad;
    plot_count = 0;
    fig_count = fig_start;
    left_margin = 0.1;  % left margin for x and y label
    bottom_margin = 0.1;
    xgap = 0.004;  % between subplots
    ygap = 0.06; 
    for cdx=1:length(sound_struct)
        if mod(cdx, plot_col*plot_row)==1
            if cdx>1
                fn_fig = fullfile(save_folder, strcat(expID, '_', suffix, '_', num2str(fig_count), '.fig'));
                savefig(fn_fig);
                fig_count = fig_count + 1;
            end
            figure(fig_count);
            set(gcf,'Position',[10 10 1200 600]);
%             tiledlayout(plot_row, plot_col, 'TileSpacing','none');
            plot_count = 0;
        end
        plot_count = plot_count + 1;
        temp = sound_struct(cdx);
        ax = subplot(plot_row, plot_col, plot_count);
%         nexttile;
        ax = gca;
        showAudioSpectrogramZZ(temp.signalNorm, temp.fs, gca(), [500,7500]);
%         showAudioSpectrogramZZ_vector(temp.signalNorm, temp.fs, gca(), [500,7500]);
%         t = datestr(temp.onsetTimeAbs, 'HH:MM:SS');
        timeInfusion = timeofday(datetime('11_00_00', 'Format','HH_mm_ss'));
        t1 = timeofday(temp.onsetTimeAbs);
        t = floor(minutes(t1-timeInfusion));
        title(ax, sprintf('t~%d', t), 'FontSize', 8);
%         timeInfusion = datetime('20221127_11_00_00', 'Format','yyyyMMdd_HH_mm_ss');
%         t1 = temp.onsetTimeAbs;
%         t = hours(t1-timeInfusion);
%         title(ax, sprintf('t~%.1f', t), 'FontSize', 8);
    %     title(num2str(cdx));
        x_tick = [0, 1];
        xticks(x_tick);
        xticklabels(round(x_tick*(temp.offsetTimeOfFile-temp.onsetTimeOfFile), 4));
        yticks(1000:2000:7000);
        yticklabels([1:2:7]);
        if mod(plot_count, plot_col)~=1
            set(ax, 'YTickLabel',[]);
        end
        set(ax,'TickDir','out');
        % change subplot location
        plot_i = ceil(plot_count/plot_col);
        plot_j = plot_count - plot_col*(plot_i-1);
%         disp([plot_i plot_j]);
        plot_x = left_margin + width*(plot_j-1) + xgap*(plot_j-1);
        plot_y = 1 - bottom_margin - height*plot_i - ygap*plot_i;
%         disp([plot_x plot_y width height]);
        set(ax, 'Position', [plot_x plot_y width height]);
        if plot_i~=plot_row
            set(ax, 'XTickLabel',[]);
        end
    end
    fn_fig = fullfile(save_folder, strcat(expID, '_', suffix, '_', num2str(fig_count), '.fig'));
    savefig(fn_fig);
end

