%Imports
import IBWread.*
import readIBWheaders.*
import readIBWbinheader.*

%Data path
datapath = '';
foldername = '';

%Experiment-specific settings
exclude_files = ["test"];
probe_type = "pyramid"; %Choose contact model to use: colloid, pyramid, or blunted pyramid

%Load data
files = struct2table(dir([datapath,foldername,'\*.ibw']));
nfiles = length(files.name);
ntiles = ceil(sqrt(nfiles));
[kc, rate] = get_parms([files.folder{2},'\',files.name{2}]);

%Constants
dwell = 2; %s
timepts = [dwell, 12,22,26.2913,29.2913,31.2913,31.7913,31.9913,32.1913];%programmed durations of each frequency, s
freqs = [0.1,0.3,0.7,1,3,10,30,50]; %frequencies tested, Hz

nu = 0.5; %assume incompressible
b0 = 5.2e-6; %nN s / nm, hydrodynamic drag coefficient
half_angle = 35*pi/180; %rad, half angle for pyramidal probe
R_coll = 12.5; %um, radius of sphere for colloidal probe
R_blunt = 40; %nm, tip radius of curvature for blunted pyramid probe

i_approach = 5000;
i_sampling = 1000;
i_times = floor(rate * timepts);
cp_tolerance = 2; %tuning parameter for detecting contact point

%Initialize results and plots
results = table();
results.FileNames = files.name;
results.Eprimes = nan(nfiles,8);
results.Edoubleprimes = nan(nfiles,8);
results.IndDepth = nan(nfiles,1);
fig1 = figure(1);
layout1 = tiledlayout(ntiles,ntiles,'TileSpacing','none','Padding','tight');
ylabel(layout1,'Raw deflection (nm)'); xlabel(layout1,'Datapoint');
fig2 = figure(2);

for I = 1:nfiles
    %Read file type
    filename = files.name{I};
    if contains(filename,exclude_files)
        continue;
    end

    %Open and plot raw data
    rawdata = IBWread([files.folder{I},'\',filename]).y;
    Z_raw = rawdata(:,1); defl_raw = rawdata(:,2); %m (Z piezo, cantilever deflection)
    figure(1); nexttile(I); hold on; plot(defl_raw);
    if length(Z_raw) < i_times(end)
        figure(1); nexttile(I); text(0.1,0.2,'Incomplete curve','Units','normalized');
        continue;
    end

    %Find contact point from initial approach
    defl_init = defl_raw(1:i_approach) * 1e9; [defl_max,i_start] = max(defl_init);
    i_relaxed = i_start + rate * dwell;
    defl_init = defl_raw(1:i_relaxed) * 1e9;
    fit_init = polyfit(1:i_sampling,defl_init(1:i_sampling),1);
    defl_init = defl_init - fit_init(1) * [1:i_relaxed]' - fit_init(2);
    threshold = mean(defl_init(1:i_sampling)) + cp_tolerance*std(defl_init(1:i_sampling));
    i_cp = i_start;
    while any(defl_init(i_cp - 1:i_cp) > threshold)
        i_cp = i_cp - 1;
    end
    if i_cp < i_sampling
        defl_init = defl_raw(1:i_relaxed) * 1e9;
        fit_init = polyfit(1:i_sampling/2,defl_init(1:i_sampling/2),1);
        defl_init = defl_init - fit_init(1) * [1:i_relaxed]' - fit_init(2);
        threshold = mean(defl_init(1:i_sampling/2)) + cp_tolerance*std(defl_init(1:i_sampling/2));
        i_cp = i_start;
        while any(defl_init(i_cp - 2:i_cp) > threshold)
            i_cp = i_cp - 1;
        end
        if i_cp < i_sampling/2
            figure(1); nexttile(I);
            text(0.1,0.2,'Insufficient baseline','Units','normalized');
            continue;
        end
    end

    Z_cp = Z_raw(i_cp) * 1e9; %contact point Z position, nm

    %Plot contact and start points
    figure(1); nexttile(I);
    plot(i_cp,defl_raw(i_cp),'ko');
    plot(i_start,defl_raw(i_start),'ro');
    text(0.1,0.1,filename,'Interpreter','none','Units','normalized');

    %Calculate geometry parameters
    Z_indented = Z_raw(i_relaxed) * 1e9 - defl_init(end); %Z position of tip after indentation & relaxation, nm
    h_indent = Z_indented - Z_cp; %total depth of indentation, nm

    %Initialize results and plots for Fourier transform analysis
    Z_fft = nan(length(freqs),3);
    F_fft = nan(length(freqs),3);
    Z_raw = Z_raw(i_start:i_start + i_times(end)) * 1e9; %Z, nm
    F_raw = defl_raw(i_start:i_start + i_times(end)) * 1e9 * kc; %force, nN

    figure(2); layout2 = tiledlayout(4,4,'TileSpacing','tight','Padding','tight');
    ylabel(layout2,'Tip position (nm) or force (nN)'); xlabel(layout2,'Time (s)');

    %Fourier transform analysis for each frequency tested
    for j = 1:length(freqs)
        time = transpose(timepts(j) : 1/rate : timepts(j + 1)) - timepts(j); %from start of frequency, s
        L = length(time);
        i_t1 = i_times(j); i_t2 = i_times(j) + L - 1;

        %Background subtract F and Z
        fit_F = polyfit(time,F_raw(i_t1:i_t2),1);
        F = detrend(F_raw(i_t1:i_t2),0);
        Z_tip = Z_raw(i_t1:i_t2) - mean(Z_raw(i_t1:i_t2)) - F / kc;

        %Zero pad to give frequency resolution of  < 0.002 Hz
        Z_pad = [Z_tip; zeros(2^20 - L,1)];
        F_pad = [F; zeros(2^20 - L,1)];
        L_pad = length(Z_pad);

        %Perform Fourier transforms, get oscillation parameters
        YZ = fft(Z_pad);
        i_approx = ceil(freqs(j) * L_pad / rate);
        [peakZ,i_freq_Z] = max(abs(YZ(i_approx - 1:i_approx + 1)));
        i_freq_Z = i_freq_Z + i_approx - 1;
        f_Z = (i_freq_Z - 1) * rate / L_pad; %fitted frequency
        amp_Z = abs(YZ(i_freq_Z) / L) * 2; %fitted amplitude
        phase_Z = angle(YZ(i_freq_Z)); %fitted phase
        Z_fft(j,:) = [amp_Z,f_Z,phase_Z]; %nm, Hz, -

        YF = fft(F_pad);
        [peakF,i_freq_F] = max(abs(YF(i_approx - 1:i_approx + 1)));
        i_freq_F = i_freq_F + i_approx - 1;
        f_F = (i_freq_F - 1) * rate / L_pad;
        amp_F = abs(YF(i_freq_F) / L) * 2;
        phase_F = angle(YF(i_freq_F));
        F_fft(j,:) = [amp_F,f_F,phase_F];

        %Plot fitting results
        nexttile(j); hold on;
        plot(time,Z_tip,'k');
        plot(time,amp_Z * cos(f_Z * time * 2 * pi + phase_Z),'LineWidth',1);
        nexttile(j + 8); hold on;
        plot(time,F,'k');
        plot(time,amp_F * cos(f_F * time * 2 * pi + phase_F),'LineWidth',1);
    end
    %Calculate compressive moduli
    switch probe_type
        case "colloid" %Mahaffy et al., 2000:
            Eprime = 1e6*(1-nu^2)/(2*sqrt(R*h_indent/1000)) * F_fft(:,1)./Z_fft(:,1) .* cos(F_fft(:,3)-Z_fft(:,3)-0);
            Eprimeprime = 1e6*(1-nu^2)/(2*sqrt(R*h_indent/1000)) * (F_fft(:,1)./Z_fft(:,1) .* sin(F_fft(:,3)-Z_fft(:,3)-0) - Z_fft(:,2)*b0);
        case "pyramid"
            Eprime = 1e9*(1-nu^2)/(sqrt(2)*tan(half_angle)*h_indent)*F_fft(:,1)./Z_fft(:,1) .* cos(F_fft(:,3)-Z_fft(:,3)-0);
            Eprimeprime = 1e9*(1-nu^2)/(sqrt(2)*tan(half_angle)*h_indent)*(F_fft(:,1)./Z_fft(:,1) .* sin(F_fft(:,3)-Z_fft(:,3)-0) - Z_fft(:,2)*b0);
        case "blunted pyramid" %Based on contact mechanics model from Rico et al., 2005:
            b = Rc*cos(half_angle); %nm, size of blunted tip region
            a = ceil(b):1:10000; %nm, effective contact radius
            d = a/tan(half_angle) * 2^(3/2) / pi .* (pi/2 - asin(b./a)) - a./Rc.*(sqrt(a.^2-b.^2) - a); %nm, indentation depth delta
            Fstar = 1e-6 * 2/(1 - nu^2) * (d.*a - sqrt(2)/pi*a.^2/tan(half_angle).*(pi/2 - asin(b./a)) - a.^3/(3*Rc) ...
                + sqrt(a.^2-b.^2).*(sqrt(2)/pi*b/tan(half_angle) + (a.^2-b.^2)/(3*Rc))); %nN/kPa. (Fstar = F/E, expression for force divided by Young's modulus)
            dFstar_da = diff(Fstar)/1; %(nN/kPa)/nm (step size of vector a is 1 nm). numerical derivative d(F/E)/da
            da_dd = 1./diff(d); %nm/nm (step size of vector a is 1 nm) numerical derivative da/d(delta)
            i_d0 = find(d>h_indent,1); %index for experimental baseline indentation depth
            a0 = a(i_d0); %experimental effective contact radius
            Eprime = 1e3 * 1/(dFstar_da(i_d0)*da_dd(i_d0)) * F_fft(:,1)./Z_fft(:,1) .* cos(F_fft(:,3)-Z_fft(:,3)-0);
            Eprimeprime = 1e3 * 1/(dFstar_da(i_d0)*da_dd(i_d0)) * (F_fft(:,1)./Z_fft(:,1) .* sin(F_fft(:,3)-Z_fft(:,3)-0) - Z_fft(:,2)*b0);
    end
    
    %Save results and figure
    results.Eprimes(I,:) = Eprime'; %Pa
    results.Edoubleprimes(I,:) = Eprimeprime'; %Pa
    results.IndDepth(I) = h_indent / 1000; %um
    saveas(figure(2),[datapath,foldername,'\',filename,'.png'])
    clf(2);
end

%save figure and export results
fig1.WindowState = 'maximized';
saveas(figure(1),[datapath,foldername,'\raw_deflections.png'])

writetable(results,'microrheology_results.csv');

function [k_c,rate] = get_parms(path)
metadata = splitlines(IBWread(path).WaveNotes);
metadata = cellfun(@(x) regexp(x,':','split','once'),metadata(1:end - 1),'UniformOutput',false);
metadata = vertcat(metadata{:});
k_c = str2double(metadata{strcmp(metadata,'SpringConstant'),2});
rate = str2double(metadata{strcmp(metadata,'NumPtsPerSec'),2});
end
