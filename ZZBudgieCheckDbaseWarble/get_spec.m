function [spec, dt, f ,T]=get_spec(signal,fs)
%     SegmentationParams;  
    params.fs=20000;
    params.min_dur=100e-3;
    params.max_dur=300e-3;
    params.nperseg=256;
    params.noverlap=128;
    params.min_fs=300;
    params.max_fs=12000;
    params.N=128;
    params.min_interval=25e-3;
    params.mel_choice=1;
    params.time_stretch=0;
    %fs=params.fs;
    nperseg=params.nperseg;
    noverlap=params.noverlap;
    min_fs=params.min_fs;
    max_fs=params.max_fs;
    min_dur=params.min_dur;
    max_dur=params.max_dur;
    N=params.N;
    mel_choice=params.mel_choice;
    time_stretch=params.time_stretch;
    interpolate_choice=0;
    % S(f,t)
    [S,F,T]=stft(signal,fs,'Window',hamming(nperseg,'periodic'),'Overlaplength',noverlap);
    ind_f=find(F>min_fs&F<max_fs);
    f=F(ind_f);
    spec=log(abs(S));
    spec=spec(ind_f,:);
    
    if interpolate_choice==1 
    %interpolation         
        if time_stretch==1
            duration = max(T)-min(T);
            duration = sqrt(duration * max_dur);
            t_stretch_interp = linspace(0, duration, length(T));
            [T_m,f_m]=meshgrid(t_stretch_interp,f);
            t_interp=linspace(0,300e-3,N);
        else
            [T_m,f_m]=meshgrid(T,f);
            t_interp=linspace(0,300e-3,N);
        end
        if mel_choice==1
            f_interp=linspace(mel(min_fs),mel(max_fs),N);
            f_interp=inv_mel(f_interp);
        else
            f_interp=linspace(min_fs,max_fs,N);
        end
        [T_interp, F_interp]=meshgrid(t_interp,f_interp);
        spec_interp=interp2(T_m,f_m,spec,T_interp,F_interp);
        % shift to have zeroes symmetrically padded
        if time_stretch==1
            indtmax=find(abs(duration-t_interp)==min(abs(duration-t_interp)));
        else
            indtmax=find(abs(max(T)-t_interp)==min(abs(max(T)-t_interp)));
        end
        shift=floor(0.5*(N-indtmax));
        spec_interp=spec_interp(:,[indtmax+shift+1:N, 1:indtmax+shift]);
        spec_interp=(spec_interp-min(min(spec_interp)))/(max(max(spec_interp))-spec_min_val);
        spec=spec_interp;
    end
    % clip value
%     spec=(spec-spec_min_val)/(spec_max_val-spec_min_val);
%     spec=clip(spec,0,1);
%     end
%     % debug
    % imagesc(t_interp,f_interp,spec_interp)
    % set(gca,'YDir','normal')
    dt=T(2)-T(1);
end
