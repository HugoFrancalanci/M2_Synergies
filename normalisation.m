clc; clear all;

% Ajout de la bibliothèque btk
addpath(genpath('C:\Users\Francalanci Hugo\Documents\MATLAB\Stage Sainte-Justine\HUG\btk'));

% Définition des sujets
subjects = {'TF02'};
nb_subjects = length(subjects);

% Définition des muscles 
muscles = {'RDELTA', 'RDELTM', 'RDELTP', 'RTRAPM', 'RTRAPS', 'RLATD', 'RSERRA'};
nb_muscles = length(muscles);

% Définition des labels pour les tâches
functional_labels = {'Functional Task 1', 'Functional Task 2', 'Functional Task 3', 'Functional Task 4'};
analytic_labels = {'Analytic Task 1', 'Analytic Task 2', 'Analytic Task 3', 'Analytic Task 4'};

nb_functional = length(functional_labels);
nb_analytic = length(analytic_labels);

% Paramètres EMG
fs = 2000;
[b, a] = butter(4, [15, 475] / (fs/2), 'bandpass');
rms_window = round(0.250 * fs);
num_points = 1000;
time_normalized = linspace(0, 1, num_points);

% Initialisation des MVCs
mvc_flexion = zeros(nb_muscles, nb_analytic);
all_functional_data = cell(nb_subjects, nb_functional, nb_muscles);

% Sélection automatique des tâches analytiques pour chaque muscle
for subj_idx = 1:nb_subjects
    for m = 1:nb_muscles
        max_activity = 0;
        best_analytic = 1;
        
        for analytic_idx = 1:nb_analytic
            fileName_analytic = sprintf(['C:\\Users\\Francalanci Hugo\\Documents\\MATLAB\\Stage Sainte-Justine\\HUG\\Sujets\\%s\\' ...
                '%s-%s-20240101-PROTOCOL01-ANALYTIC%d-.c3d'], ...
                subjects{subj_idx}, subjects{subj_idx}, subjects{subj_idx}, analytic_idx);
            
            c3dH_analytic = btkReadAcquisition(fileName_analytic);
            analogs_analytic = btkGetAnalogs(c3dH_analytic);
            
            muscle_name = muscles{m};
            if isfield(analogs_analytic, muscle_name)
                signal = analogs_analytic.(muscle_name);
                signal_filtered = filtfilt(b, a, signal);
                signal_abs = abs(signal_filtered);
                emg_rms = sqrt(movmean(signal_abs.^2, rms_window));
                mvc_flexion(m, analytic_idx) = max(emg_rms);
                
                if mvc_flexion(m, analytic_idx) > max_activity
                    max_activity = mvc_flexion(m, analytic_idx);
                    best_analytic = analytic_idx;
                end
            end
        end
        selected_analytics(m) = best_analytic;
    end
end

% Sélection d'une tâche fonctionnelle unique
selected_functional = input('Choisissez une tâche fonctionnelle (1-4): ');

% Traitement et normalisation de la tâche fonctionnelle choisie
for subj_idx = 1:nb_subjects
    fileName_functional = sprintf(['C:\\Users\\Francalanci Hugo\\Documents\\MATLAB\\Stage Sainte-Justine\\HUG\\Sujets\\%s\\' ...
        '%s-%s-20240101-PROTOCOL01-FUNCTIONAL%d-.c3d'], ...
        subjects{subj_idx}, subjects{subj_idx}, subjects{subj_idx}, selected_functional);
    
    c3dH_functional = btkReadAcquisition(fileName_functional);
    analogs_functional = btkGetAnalogs(c3dH_functional);
    
    figure;
    sgtitle(sprintf('Sujet %s - %s', subjects{subj_idx}, functional_labels{selected_functional}));
    
    for m = 1:nb_muscles
        subplot(ceil(nb_muscles/2), 2, m);
        muscle_name = muscles{m};
        if isfield(analogs_functional, muscle_name)
            signal = analogs_functional.(muscle_name);
            signal_filtered = filtfilt(b, a, signal);
            signal_abs = abs(signal_filtered);
            emg_rms = sqrt(movmean(signal_abs.^2, rms_window));
            emg_normalized = (emg_rms / mvc_flexion(m, selected_analytics(m))) * 100;
            interp_signal_functional = interp1(linspace(0, 1, length(signal)), emg_normalized, time_normalized, 'spline');
            all_functional_data{subj_idx, selected_functional, m} = interp_signal_functional;
            plot(time_normalized, interp_signal_functional, 'LineWidth', 1.2);
            title(muscle_name);
            ylabel('EMG (% MVC)');
            grid on;
        else
            title(sprintf('%s (Données absentes)', muscle_name));
        end
    end
end

