clear
%% Загрузка исходных данных и настроек
run('settings_RTDS_runtime_PMU_RC');

%% Расчет недостающих данных
% Параметры ВЛ
HVL_params.z0 = HVL_params.r0 + HVL_params.x0 * 1i; % Погонное полное сопротивление нулевой последовательности ВЛ
HVL_params.z1 = HVL_params.r1 + HVL_params.x1 * 1i; % Погонное полное сопротивление прямой последовательности ВЛ

% Параметры расчета
if calc_settings.need_factors_ranges == true
    
    % Нагрузка
    calc_settings.load.exp_num = round((calc_settings.load.finish_value - ...
        calc_settings.load.start_value) / calc_settings.load.step_change) + 1;
    calc_settings.load.range = calc_settings.load.start_value:...
        calc_settings.load.step_change:calc_settings.load.finish_value;

    % Местоположение КЗ
    calc_settings.sc_position.exp_num = round((calc_settings.sc_position.finish_value - ...
        calc_settings.sc_position.start_value) / calc_settings.sc_position.step_change) + 1;
    calc_settings.sc_position.range = calc_settings.sc_position.start_value:...
        calc_settings.sc_position.step_change:calc_settings.sc_position.finish_value;

    % Фаза КЗ
    calc_settings.sc_phase.exp_num = round((calc_settings.sc_phase.finish_value - ...
        calc_settings.sc_phase.start_value) / calc_settings.sc_phase.step_change) + 1;
    calc_settings.sc_phase.range = calc_settings.sc_phase.start_value:...
        calc_settings.sc_phase.step_change:calc_settings.sc_phase.finish_value;

    % Взаимный угол ЭДС
    calc_settings.ef_mutual_angle.exp_num = round((calc_settings.ef_mutual_angle.finish_value - ...
        calc_settings.ef_mutual_angle.start_value) / calc_settings.ef_mutual_angle.step_change) + 1;
    calc_settings.ef_mutual_angle.range = calc_settings.ef_mutual_angle.start_value:...
        calc_settings.ef_mutual_angle.step_change:calc_settings.ef_mutual_angle.finish_value;

    % Переходное сопротивление
    calc_settings.trans_resistance.exp_num = ...
        round((calc_settings.trans_resistance.finish_value_int1 - ...
        calc_settings.trans_resistance.start_value_int1) / ...
        calc_settings.trans_resistance.step_change_int1) + 1 + ...
        round((calc_settings.trans_resistance.finish_value_int2 - ...
        calc_settings.trans_resistance.start_value_int2) / ...
        calc_settings.trans_resistance.step_change_int2) + 1 + ...
        round((calc_settings.trans_resistance.finish_value_int3 - ...
        calc_settings.trans_resistance.start_value_int3) / ...
        calc_settings.trans_resistance.step_change_int3) + 1;
    calc_settings.trans_resistance.range = [...
        calc_settings.trans_resistance.start_value_int1:...
        calc_settings.trans_resistance.step_change_int1:...
        calc_settings.trans_resistance.finish_value_int1 ...
        calc_settings.trans_resistance.start_value_int2:...
        calc_settings.trans_resistance.step_change_int2:...
        calc_settings.trans_resistance.finish_value_int2 ...
        calc_settings.trans_resistance.start_value_int3:...
        calc_settings.trans_resistance.step_change_int3:...
        calc_settings.trans_resistance.finish_value_int3 ...
        ];

    % Отношение сопротивления ПП системы справа к системе слева
    calc_settings.positive_seq.exp_num = round((calc_settings.positive_seq.finish_value - ...
        calc_settings.positive_seq.start_value) / calc_settings.positive_seq.step_change) + 1;
    calc_settings.positive_seq.range = calc_settings.positive_seq.start_value:...
        calc_settings.positive_seq.step_change:calc_settings.positive_seq.finish_value;

    % Отношение сопротивления НП системы справа к системе слева
    calc_settings.zero_seq.exp_num = round((calc_settings.zero_seq.finish_value - ...
        calc_settings.zero_seq.start_value) / calc_settings.zero_seq.step_change) + 1;
    calc_settings.zero_seq.range = calc_settings.zero_seq.start_value:...
        calc_settings.zero_seq.step_change:calc_settings.zero_seq.finish_value;

end

%% Основная прграмма
for factor = calc_settings.factors

    %% Загрузка данных УСВИ
    if calc_settings.need_load_PMU_from_csv == true
        if PMU_settings.is_not_runtime_PMU == true
            PMU_data_matrix = load_PMU_data_from_csv(...
                directories_settings, factor, calc_settings, ...
                PMU_settings, SV_settings, sc_settings);
        else     
            PMU_data_matrix = load_RTDS_runtime_PMU_data_from_csv(...
                directories_settings, factor, calc_settings);
        end
    else
        init_data_file_name = strcat('PMU_data_matrix_',factor);
        init_data_path = fullfile(directories_settings.PMU_data_MAT, init_data_file_name);
        load(init_data_path);
    end

    %% Получение данный о токах и напряжениях с преобразование в комплексные
    % числа из данных PMU
    [Ibeg, Iend, Ubeg, Uend] = get_IU_from_PMU(PMU_data_matrix, factor, ...
        PMU_settings, calc_settings);    
    
    %% Вычисление погрешностей выбранного алгоритма во временной области
    [delta_percent_time_zone] = calc_algs_error_time_zone(...
        Ibeg, Iend, Ubeg, Uend, factor, HVL_params, ...
        calc_settings, sc_settings, directories_settings);
    
    %% Определение моментов КЗ для случая изменения фазы КЗ
    if factor == "sc_phase"
        sc_settings.fault_inception_moment.sc_phase = ...
            fault_inception_determening(...
            factor, calc_settings, PMU_settings, sc_settings);
    end
    
    %% Построение зависимостей погрешностей от комплекта PMU
    create_errors_graphs_time_zone(...
        delta_percent_time_zone, factor, ...
        PMU_settings.period, sc_settings, graphs_settings, calc_settings);
    
    clear delta_percent;
    
    %% Вычисление погрешностей алгоритмов для конкретного комплекта PMU
    delta_percent = calc_algs_error_one_set_PMU(...
        Ibeg, Iend, Ubeg, Uend, factor, PMU_settings.period, HVL_params, ...
        calc_settings, sc_settings, directories_settings);

    %% Построение зависимостей погрешностей для конкретного комплекта PMU
    create_errors_graphs_one_set_PMU(...
        delta_percent, factor, calc_settings, graphs_settings);
    
    clear delta_percent;
    
    if factor == "sc_phase"
        delta_percent = Monte_Carlo_line_params_variation(...
            Ibeg, Iend, Ubeg, Uend, factor, PMU_settings.period, HVL_params, ...
            calc_settings, sc_settings, directories_settings);
    end
    
    %% Освобождение памяти от лишних переменных
    clearvars -except HVL_params calc_settings sc_settings PMU_settings ...
        graphs_settings graphs_settings directories_settings;
    
end

clear;

%% Функции

% Функция загрузки данных runtime УСВИ RTDS из CSV
function PMU_data_matrix = load_RTDS_runtime_PMU_data_from_csv(directories_settings, factor, calc_settings)
    exp_count = getfield(calc_settings,factor).exp_num;
    
    for idx = 1:exp_count
        file_name = strcat("PMU_data_", factor, "_exp", num2str(idx), ".csv");
        file_path = fullfile(directories_settings.PMU_data_CSV, file_name);
        PMU_data_cells{idx} = csvread(file_path,1,0);
    end
    
    PMU_data_matrix = cell2mat(PMU_data_cells);
    
    results_file_name = strcat("PMU_data_matrix_",factor, ".mat");
    results_path = fullfile(directories_settings.PMU_data_MAT, results_file_name);
    save(results_path, 'PMU_data_matrix');

end

% Функция загрузки данных УСВИ из CSV
function PMU_data_matrix = load_PMU_data_from_csv(...
    directories_settings, factor, calc_settings, PMU_settings, SV_settings, sc_settings)

    file_name = strcat("PMU_data_beg_line_", factor, ".csv");
    file_path = fullfile(directories_settings.PMU_data_CSV, file_name);
    names = getfield(PMU_settings.names,PMU_settings.PMU_data_beg_line);
    PMU_data_beg_line = convert_csv_to_timetable(file_path,names);
    
    file_name = strcat("PMU_data_end_line_", factor, ".csv");
    file_path = fullfile(directories_settings.PMU_data_CSV, file_name);
    names = getfield(PMU_settings.names,PMU_settings.PMU_data_end_line);
    PMU_data_end_line = convert_csv_to_timetable(file_path,names);   

    file_name = strcat("SV_beg_line_", factor, ".csv");
    file_path = fullfile(directories_settings.SV_data_CSV, file_name);
    names = SV_settings.names.SV_beg_line;
    SV_beg_line = convert_csv_to_timetable(file_path,names);
    SV_beg_line = parse_SV_256(SV_beg_line,names);
    
    file_name = strcat("SV_end_line_", factor, ".csv");
    file_path = fullfile(directories_settings.SV_data_CSV, file_name);
    names = SV_settings.names.SV_end_line;
    SV_end_line = convert_csv_to_timetable(file_path,names);
    SV_end_line = parse_SV_256(SV_end_line,names);

    plot(...
        PMU_data_beg_line.Timestamp,PMU_data_beg_line.IphsAmag,...
        PMU_data_end_line.Timestamp,PMU_data_end_line.IphsAmag,...
        SV_beg_line.Timestamp-SV_settings.delta_time,SV_beg_line.IphsA,...
        SV_end_line.Timestamp-SV_settings.delta_time,SV_end_line.IphsA);
    
    PMU_data_matrix = identify_short_circuit(...
        SV_beg_line, SV_end_line, PMU_data_beg_line, PMU_data_end_line, ...
        calc_settings, SV_settings, PMU_settings, sc_settings); 
    
%     for idx = 1:exp_count
%         file_name = strcat("PMU_data_", factor, "_exp", num2str(idx), ".csv");
%         file_path = fullfile(directories_settings.PMU_data_CSV, file_name);
%         PMU_data_cells{idx} = csvread(file_path,1,0);
%     end
%     
%     PMU_data_matrix = cell2mat(PMU_data_cells);
    
    results_file_name = strcat("PMU_data_matrix_",factor, ".mat");
    results_path = fullfile(directories_settings.PMU_data_MAT, results_file_name);
    save(results_path, 'PMU_data_matrix');

end

% Функция чтения csv и записи его данных в timetable
function table = convert_csv_to_timetable(file_path, names)
    opts = detectImportOptions(file_path, 'Delimiter', ',');
    table = readtimetable(file_path, opts);
    table.Properties.DimensionNames{1} = 'Timestamp';
    for idx=1:length(names)
        table.Properties.VariableNames{idx} = names{1,idx};
    end
end

% Функция разбора данных для SV с 256 выборками за период пром. частоты
function OutTable = parse_SV_256(InTable,names)
    multiplier_current = 0.001;
    multiplier_voltage = 0.01;
    T = '00:00.000078125';
    infmt = 'mm:ss.SSSSSSSSS';
    dt = duration(T,'InputFormat',infmt);
    values = zeros(length(InTable.Timestamp)*8,8);
    time = NaT(length(InTable.Timestamp)*8,1);
    for idx_packet=1:8
        time(idx_packet:8:end) = InTable.Timestamp  - dt * (8 - idx_packet);
        for idx_value=1:4
            values(idx_packet:8:end,idx_value) = InTable(:,idx_value+8*(idx_packet-1)).Variables*multiplier_current;
        end
        for idx_value=5:8
            values(idx_packet:8:end,idx_value) = InTable(:,idx_value+8*(idx_packet-1)).Variables*multiplier_voltage;
        end
    end
    OutTable = timetable(values(:,1),values(:,2),values(:,3),values(:,4),values(:,5),values(:,6),values(:,7),values(:,8),'RowTimes',time);
    OutTable.Properties.DimensionNames{1} = 'Timestamp';
    for idx=1:length(names)
        OutTable.Properties.VariableNames{idx} = names{1,idx};
    end
end

% Функция определения КЗ на осциллограмме
function PMU_data_matrix = identify_short_circuit(...
    SV_beg_line, SV_end_line, PMU_data_beg_line, PMU_data_end_line, ...
    calc_settings, SV_settings, PMU_settings, sc_settings)

%     idx_current_beg=find(abs(SV_beg_line(:,1).Variables)>=calc_settings.threshold_val,1,'first');
%     while idx_current_beg > 2
%         if (SV_beg_line(idx_current_beg,1).Variables >= 0 && ...
%                 SV_beg_line(idx_current_beg-1,1).Variables <= 0) ...
%                 || (SV_beg_line(idx_current_beg,1).Variables <= 0 ...
%                 && SV_beg_line(idx_current_beg-1,1).Variables >= 0)
%             current_beg_time = SV_beg_line.Timestamp(idx_current_beg) - ...
%                 SV_settings.delta_time;
%             break;
%         end 
%         idx_current_beg = idx_current_beg - 1;
%     end

%     sc_beg_time = current_beg_time + sc_settings.prefault_duration;
%     sc_end_time = sc_beg_time + sc_settings.fault_duration;
% 
%     idx_sc_beg=find(PMU_data_beg_line.Timestamp>=(sc_beg_time - PMU_settings.PMU_shift_beg_line),1,'first');
%     idx_sc_end=find(PMU_data_beg_line.Timestamp>=(sc_end_time - PMU_settings.PMU_shift_end_line),1,'first');

    sc_num = 0;
    idx_current_beg = find(PMU_data_beg_line.IphsAmag>=calc_settings.threshold_val,1,'first');
    while (idx_current_beg)
        sc_num = sc_num + 1;
        while idx_current_beg > 1
            if (PMU_data_beg_line.IphsAmag(idx_current_beg,1) >= 0 && ...
                    PMU_data_beg_line.IphsAmag(idx_current_beg-1,1) <= 0) ...
                    || (PMU_data_beg_line.IphsAmag(idx_current_beg,1) <= 0 ...
                    && PMU_data_beg_line.IphsAmag(idx_current_beg-1,1) >= 0)
                current_beg_time = PMU_data_beg_line.Timestamp(idx_current_beg);
                break;
            end 
            idx_current_beg = idx_current_beg - 1;
        end

        sc_beg_time = current_beg_time + sc_settings.prefault_duration;
        sc_end_time = sc_beg_time + sc_settings.fault_duration;

        idx_sc_beg_PMU_data_beg_line=find(...
            PMU_data_beg_line.Timestamp>=(sc_beg_time+PMU_settings.PMU_shift_beg_line),1,'first');
        idx_sc_end_PMU_data_beg_line=find(...
            PMU_data_beg_line.Timestamp>=(sc_end_time+PMU_settings.PMU_shift_beg_line),1,'first');

        idx_sc_beg_PMU_data_end_line=find(...
            PMU_data_end_line.Timestamp>=(sc_beg_time+PMU_settings.PMU_shift_beg_line),1,'first');
        idx_sc_end_PMU_data_end_line=find(...
            PMU_data_end_line.Timestamp>=(sc_end_time+PMU_settings.PMU_shift_beg_line),1,'first');

        PMU_data_matrix(:,1+(sc_num-1)*24) = ...
            PMU_data_beg_line.IphsAmag(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,2+(sc_num-1)*24) = ...
            PMU_data_beg_line.IphsAan(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,3+(sc_num-1)*24) = ...
            PMU_data_beg_line.IphsBmag(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,4+(sc_num-1)*24) = ...
            PMU_data_beg_line.IphsBan(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,5+(sc_num-1)*24) = ...
            PMU_data_beg_line.IphsCmag(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,6+(sc_num-1)*24) = ...
            PMU_data_beg_line.IphsCan(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);

        PMU_data_matrix(:,7+(sc_num-1)*24) = ...
            PMU_data_end_line.IphsAmag(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,8+(sc_num-1)*24) = ...
            PMU_data_end_line.IphsAan(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,9+(sc_num-1)*24) = ...
            PMU_data_end_line.IphsBmag(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,10+(sc_num-1)*24) = ...
            PMU_data_end_line.IphsBan(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,11+(sc_num-1)*24) = ...
            PMU_data_end_line.IphsCmag(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,12+(sc_num-1)*24) = ...
            PMU_data_end_line.IphsCan(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);

        PMU_data_matrix(:,13+(sc_num-1)*24) = ...
            PMU_data_beg_line.VphsAmag(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,14+(sc_num-1)*24) = ...
            PMU_data_beg_line.VphsAan(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,15+(sc_num-1)*24) = ...
            PMU_data_beg_line.VphsBmag(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,16+(sc_num-1)*24) = ...
            PMU_data_beg_line.VphsBan(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,17+(sc_num-1)*24) = ...
            PMU_data_beg_line.VphsCmag(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);
        PMU_data_matrix(:,18+(sc_num-1)*24) = ...
            PMU_data_beg_line.VphsCan(idx_sc_beg_PMU_data_beg_line:idx_sc_end_PMU_data_beg_line,1);

        PMU_data_matrix(:,19+(sc_num-1)*24) = ...
            PMU_data_end_line.VphsAmag(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,20+(sc_num-1)*24) = ...
            PMU_data_end_line.VphsAan(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,21+(sc_num-1)*24) = ...
            PMU_data_end_line.VphsBmag(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,22+(sc_num-1)*24) = ...
            PMU_data_end_line.VphsBan(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,23+(sc_num-1)*24) = ...
            PMU_data_end_line.VphsCmag(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);
        PMU_data_matrix(:,24+(sc_num-1)*24) = ...
            PMU_data_end_line.VphsCan(idx_sc_beg_PMU_data_end_line:idx_sc_end_PMU_data_end_line,1);

        idx_current_beg=find(PMU_data_beg_line.IphsAmag(idx_current_beg+...
            calc_settings.vectors_between_sc:end,1)>=...
            calc_settings.threshold_val,1,'first') + ...
            idx_current_beg+calc_settings.vectors_between_sc;
    end

end

% Функция выделения токов и напряжений из матрицы данных PMU
function [Ibeg, Iend, Ubeg, Uend] = get_IU_from_PMU(PMU, factor, PMU_settings, calc_settings)

    sign_Iend = 1;
    if PMU_settings.need_change_direction_Iend == true
        sign_Iend = -1;
    end
    
    angle_coeff = 1;
    if PMU_settings.angle_in_degrees == true
        angle_coeff = pi/180;
    end
    
    RC_angle_coef = 0;
    RC_rms_coef = 1;
    if PMU_settings.RC_data == true
        RC_angle_coef = -pi/2;
        RC_rms_coef = 1000/(2 * pi * 50);
    end

    exp_count = getfield(calc_settings,factor).exp_num;   

    for idx_exp = 1:exp_count
        Ibeg.phase_A(:,idx_exp) = RC_rms_coef * ...
            PMU(:,PMU_settings.Ibeg.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Ibeg.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count)))) * ...
            exp(1j*RC_angle_coef);
        
        Ibeg.phase_B(:,idx_exp) = RC_rms_coef * ...
            PMU(:,PMU_settings.Ibeg.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Ibeg.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count)))) * ...
            exp(1j*RC_angle_coef);
        
        Ibeg.phase_C(:,idx_exp) = RC_rms_coef * ...
            PMU(:,PMU_settings.Ibeg.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Ibeg.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count)))) * ...
            exp(1j*RC_angle_coef);

        Iend.phase_A(:,idx_exp) = RC_rms_coef * sign_Iend * ...
            PMU(:,PMU_settings.Iend.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Iend.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count)))) * ...
            exp(1j*RC_angle_coef);
        
        Iend.phase_B(:,idx_exp) = RC_rms_coef * sign_Iend * ...
            PMU(:,PMU_settings.Iend.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Iend.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count)))) * ...
            exp(1j*RC_angle_coef);
        
        Iend.phase_C(:,idx_exp) = RC_rms_coef * sign_Iend * ...
            PMU(:,PMU_settings.Iend.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Iend.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count)))) * ...
            exp(1j*RC_angle_coef);
      
        Ubeg.phase_A(:,idx_exp) = PMU(:,PMU_settings.Ubeg.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Ubeg.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Ubeg.phase_B(:,idx_exp) = PMU(:,PMU_settings.Ubeg.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Ubeg.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Ubeg.phase_C(:,idx_exp) = PMU(:,PMU_settings.Ubeg.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Ubeg.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));

        Uend.phase_A(:,idx_exp) = PMU(:,PMU_settings.Uend.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Uend.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Uend.phase_B(:,idx_exp) = PMU(:,PMU_settings.Uend.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Uend.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Uend.phase_C(:,idx_exp) = PMU(:,PMU_settings.Uend.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*angle_coeff*(PMU(:,PMU_settings.Uend.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
    end
end

% Функция определения моментов КЗ для случая варьирования фазы замыкания
function fault_inception_moments = fault_inception_determening(...
    factor, calc_settings, PMU_settings, sc_settings)
    
    exp_count = getfield(calc_settings,factor).exp_num;
    
    fault_inception_moments = zeros(exp_count, 1);
    fault_inception_moments(1,1) = getfield(sc_settings.fault_inception_moment,factor);
    period = PMU_settings.period;
    for idx = 2:exp_count
        fault_inception_moments(idx,1) = fault_inception_moments(1,1) + ...
            round(period * (idx / exp_count));
    end    
end

% Функция вычисления погрешностей алгоритма во временной зоне (для всех
% комплектов PMU)
function [delta_percent] = calc_algs_error_time_zone(...
    Ibeg, Iend, Ubeg, Uend, factor, HVL_params, ...
    calc_settings, sc_settings, directories_settings)

    exp_count = getfield(calc_settings,factor).exp_num;
    
    PMU_set_count = size(Ibeg.phase_A, 1);
    delta_percent = zeros(exp_count, PMU_set_count);
    
    sc_position = sc_settings.position_init;

    for idx_exp = 1:exp_count 
        for PMU_set = 1:PMU_set_count         
            % Формирование наборов PMU с трехфазными измерениями
            Ibeg_PMU_set = [
                Ibeg.phase_A(PMU_set,idx_exp); ...
                Ibeg.phase_B(PMU_set,idx_exp); ...
                Ibeg.phase_C(PMU_set,idx_exp)];
            
            Iend_PMU_set = [
                Iend.phase_A(PMU_set,idx_exp); ...
                Iend.phase_B(PMU_set,idx_exp); ...
                Iend.phase_C(PMU_set,idx_exp)];

            Ubeg_PMU_set = [
                Ubeg.phase_A(PMU_set,idx_exp); ...
                Ubeg.phase_B(PMU_set,idx_exp); ...
                Ubeg.phase_C(PMU_set,idx_exp)];
            
            Uend_PMU_set = [
                Uend.phase_A(PMU_set,idx_exp); ...
                Uend.phase_B(PMU_set,idx_exp); ...
                Uend.phase_C(PMU_set,idx_exp)];

            % Вычисления
            if factor == "sc_position"
                sc_position = HVL_params.length * ...
                    calc_settings.sc_position.range(idx_exp);
            end
                
            dist_km = ...
                FAULT_LOCATION_EXPRESSIONS_CORRECTED_NEW_2021 ...
                (Ubeg_PMU_set, Ibeg_PMU_set, Uend_PMU_set, Iend_PMU_set, ...
                HVL_params.z1, HVL_params.z0, HVL_params.y1, HVL_params.y0, ...
                HVL_params.length, sc_position, calc_settings.alg_num, sc_settings.faulted_phase);

                delta_percent(idx_exp, PMU_set) = ...
                100 * ((real(dist_km) - sc_position) ...
                / HVL_params.length);                
        end
    end
    
    results_file_name = strcat("delta_time_zone_",factor, ".mat");
    results_path = fullfile(directories_settings.results_data, results_file_name);
    save(results_path, 'delta_percent', 'factor', 'calc_settings');
end

% Функция построения графиков погрешностей алгоритмов во временной зоне
function create_errors_graphs_time_zone(...
    delta_percent, factor, period, sc_settings, graphs_settings, calc_settings)

    exp_count = getfield(calc_settings,factor).exp_num;
    fault_inception_moment = getfield(sc_settings.fault_inception_moment,factor);

    % Построение графиков и вертикальных прямых (момент КЗ, 1 период, 2 периода)
    graph_directory = graphs_settings.result_directory;
        
    if factor == "sc_phase"
        for idx_exp = 1:exp_count
            fig = figure;
            plot(real(abs(delta_percent(idx_exp,:))))
            xlim([graphs_settings.xlim_min graphs_settings.xlim_max])
            ylim([graphs_settings.ylim_min graphs_settings.ylim_max])
            
            xlabel('Measurement number');
            ylabel('Error, %');
            
            title('Зависимость погрешности ОМП от комплекта.Изменение фазы КЗ','FontSize',8)
            
            hold on;
            xline(fault_inception_moment(idx_exp),'-',{'Fault inception'});
            hold on;
            xline(fault_inception_moment(idx_exp) + period,'-',{'1st period end'});
            hold on;
            xline(fault_inception_moment(idx_exp) + 2*period,'-',{'2nd period end'});
            hold on;
            xline(fault_inception_moment(idx_exp) + 3*period,'-',{'3rd period end'});   
            
            hold on
            plot([0 8000],[1 1])
            
            if graphs_settings.need_save == 1
                graph_file_name = strcat('sc_phase_time_zone_exp', num2str(idx_exp));
                results_path = fullfile(graph_directory, graph_file_name);
                saveas(gcf, results_path, 'png' );
                saveas(gcf, results_path, 'fig' );
%                 saveas(gcf, results_path, 'emf' );
            end
            
            close(fig);
        end
    else
        fig = figure;
        for idx_exp = 1:exp_count
            plot(real(abs(delta_percent(idx_exp,:))))
            xlim([graphs_settings.xlim_min graphs_settings.xlim_max])
            hold on
        end
        
        hold on;
        xline(fault_inception_moment,'-',{'Fault inception'});
        hold on;
        xline(fault_inception_moment + period,'-',{'1st period end'});
        hold on;
        xline(fault_inception_moment + 2*period,'-',{'2nd period end'});
        hold on;
        xline(fault_inception_moment + 3*period,'-',{'3rd period end'});
        
        hold on;
        plot([0 8000],[1 1]);
        
        % Наименование осей
        xlabel('Measurement number');
        ylabel('Error, %');
        
        % Наименование графика
        switch factor
            case "load"
                title('Зависимость погрешности ОМП от комплекта.Изменение нагрузки','FontSize',8)
            case "sc_position"
                title('Зависимость погрешности ОМП от комплекта.Изменение точки КЗ','FontSize',8) 
            case "ef_mutual_angle"
                title('Зависимость погрешности ОМП от комплекта.Изменение взаимного угла ЭДС систем','FontSize',8) 
            case "trans_resistance"
                title('Зависимость погрешности ОМП от комплекта.Изменение переходного сопротивления','FontSize',8)   
            case "positive_seq"
                title('Зависимость погрешности ОМП от комплекта.Изменение парам. обр. посл. системы','FontSize',8) 
            case "zero_seq"
                title('Зависимость погрешности ОМП от комплекта.Изменение парам. прямой посл. системы','FontSize',8)   
        end
        
        if graphs_settings.need_save == 1
            graph_file_name = strcat(factor, '_time_zone');
            results_path = fullfile(graph_directory, graph_file_name);
            saveas(gcf, results_path, 'png' );
            saveas(gcf, results_path, 'fig' );
%             saveas(gcf, results_path, 'emf' );
        end
        
        close(fig)
    end
end

% Функция вычисления погрешностей алгоритмов для одного комплекта PMU
function [delta_percent] = calc_algs_error_one_set_PMU(...
    Ibeg, Iend, Ubeg, Uend, factor, period, HVL_params, ...
    calc_settings, sc_settings, directories_settings)

    fault_inception_moment = getfield(sc_settings.fault_inception_moment,factor);

    if factor ~= "load" && factor ~= "sc_phase" && factor ~= "ef_mutual_angle"
        PMU_set = fault_inception_moment + ...
            calc_settings.calc_period_num * period;
    end
    
    sc_position = sc_settings.position_init;
    exp_count = getfield(calc_settings,factor).exp_num;
    
    delta_percent = zeros(calc_settings.alg_count, exp_count);

    for idx_exp = 1:exp_count
        
        if factor == "load" || factor == "sc_phase" || factor == "ef_mutual_angle"
            PMU_set = fault_inception_moment(idx_exp) + ...
                calc_settings.calc_period_num * period;
        end
        
        % Формирование наборов PMU с трехфазными измерениями
        Ibeg_PMU_set = [
            Ibeg.phase_A(PMU_set,idx_exp); ...
            Ibeg.phase_B(PMU_set,idx_exp); ...
            Ibeg.phase_C(PMU_set,idx_exp)];

        Iend_PMU_set = [
            Iend.phase_A(PMU_set,idx_exp); ...
            Iend.phase_B(PMU_set,idx_exp); ...
            Iend.phase_C(PMU_set,idx_exp)];

        Ubeg_PMU_set = [
            Ubeg.phase_A(PMU_set,idx_exp); ...
            Ubeg.phase_B(PMU_set,idx_exp); ...
            Ubeg.phase_C(PMU_set,idx_exp)];

        Uend_PMU_set = [
            Uend.phase_A(PMU_set,idx_exp); ...
            Uend.phase_B(PMU_set,idx_exp); ...
            Uend.phase_C(PMU_set,idx_exp)];

        % Вычисления
        if factor == "sc_position"
            sc_position = HVL_params.length * ...
                (calc_settings.sc_position.start_value + ...
                calc_settings.sc_position.step_change * ...
                (idx_exp - 1));
            sc_position = HVL_params.length * ...
                calc_settings.sc_position.range(idx_exp);
        end
        
        dist_km = ...
            FAULT_LOCATION_EXPRESSIONS_CORRECTED_NEW_2021 ...
            (Ubeg_PMU_set, Ibeg_PMU_set, Uend_PMU_set, Iend_PMU_set, ...
            HVL_params.z1, HVL_params.z0, HVL_params.y1, HVL_params.y0, ...
            HVL_params.length, sc_position, 1:calc_settings.alg_count, sc_settings.faulted_phase);
    
        for idx_alg = 1:calc_settings.alg_count
            delta_percent(idx_alg,idx_exp) = 100 * (real(dist_km(idx_alg)) - sc_position) / HVL_params.length;
        end
    end
    
    results_file_name = strcat('delta_',factor, ".mat");
    results_path = fullfile(directories_settings.results_data, results_file_name);
    save(results_path, 'delta_percent', 'factor', 'calc_settings');
end

% Построение графиков погрешностей алгоритмов для одного комплекта PMU
function create_errors_graphs_one_set_PMU(delta_percent, factor, calc_settings, graphs_settings)
    
    % Настройки графика
    marker_size = graphs_settings.marker_size;
    markers = graphs_settings.markers;

    range = getfield(calc_settings,factor).range;
    
    % Построение графика
    graph_directory = graphs_settings.result_directory;
    fig = figure; 
    idx = 1;
    for idx_alg = graphs_settings.algs_for_plot
        scatter(range, delta_percent(idx_alg,:), ...
            marker_size, 'Marker', markers(idx), 'MarkerEdgeColor','k');
        hold on;
        labels(idx) = strcat("Alg.",num2str(idx_alg));
        idx = idx + 1;
    end
    
    % Легенда графика
    legend (labels,'Location','best')
    
    % Подпись оси x
    switch factor
        case "load"
            xlabel('Load, p.u.');
        case "sc_position"
            xlabel('Distance from the left busbar to fault, p.u.');
        case "sc_phase"
            xlabel('Fault inception angle, deg.');
        case "ef_mutual_angle"
            xlabel('Phase angle difference between the left side and right side, deg.');
        case "trans_resistance"
            xlabel('Fault resistance, Ohm');   
        case "positive_seq"
            xlabel('Left-to-right system positive sequence resistance ratio, p.u.');
        case "zero_seq"
            xlabel('Left-to-right system zero sequence resistance ratio, p.u.');   
    end
    
    xticks(range);

    % Подпись оси y
    ylabel('Abs value of fault location error, %');
    
    axis padded;
    
    if graphs_settings.need_save == 1
        graph_file_name = strcat(factor, "_PMU_set");
        results_path = fullfile(graph_directory, graph_file_name);
        saveas(gcf, results_path, 'png' );
        saveas(gcf, results_path, 'fig' );
        saveas(gcf, results_path, 'emf' );
    end 
        
    close(fig);
end

% Функция расчета погрешностей ОМП по методу Монте-Карло при вариации
% параметров линии
function delta_percent = Monte_Carlo_line_params_variation(...
    Ibeg, Iend, Ubeg, Uend, factor, period, HVL_params, ...
    calc_settings, sc_settings, directories_settings)
    
    sc_position = sc_settings.position_init;
    fault_inception_moment = getfield(sc_settings.fault_inception_moment,factor);
    PMU_set = fault_inception_moment + calc_settings.calc_period_num * period;
    
    if factor ~= "load" && factor ~= "sc_phase" && factor ~= "ef_mutual_angle"
        PMU_set = fault_inception_moment + calc_settings.calc_period_num * period;
    else
        PMU_set = fault_inception_moment(1) + calc_settings.calc_period_num * period;
    end
    
    idx_exp = 1;
    
    % Формирование наборов PMU с трехфазными измерениями
    Ibeg_PMU_set = [
        Ibeg.phase_A(PMU_set,idx_exp); ...
        Ibeg.phase_B(PMU_set,idx_exp); ...
        Ibeg.phase_C(PMU_set,idx_exp)];

    Iend_PMU_set = [
        Iend.phase_A(PMU_set,idx_exp); ...
        Iend.phase_B(PMU_set,idx_exp); ...
        Iend.phase_C(PMU_set,idx_exp)];

    Ubeg_PMU_set = [
        Ubeg.phase_A(PMU_set,idx_exp); ...
        Ubeg.phase_B(PMU_set,idx_exp); ...
        Ubeg.phase_C(PMU_set,idx_exp)];

    Uend_PMU_set = [
        Uend.phase_A(PMU_set,idx_exp); ...
        Uend.phase_B(PMU_set,idx_exp); ...
        Uend.phase_C(PMU_set,idx_exp)];
    
    delta_percent = zeros(1000, calc_settings.alg_count);
    for idx_variation = 1:1000
        len = (0.98 + (1.02 - 0.98) * rand) * HVL_params.length;
        r1 = (0.7 + (1.3 - 0.7) * rand) * HVL_params.r1;
        r0 = (0.7 + (1.3 - 0.7) * rand) * HVL_params.r0;
        x1 = (0.98 + (1.02 - 0.98) * rand) * HVL_params.x1;
        x0 = (0.8 + (1.2 - 0.8) * rand) * HVL_params.x0;
        y1 = (0.94 + (1.06 - 0.94) * rand) * HVL_params.y1;
        y0 = (0.9 + (1.1 - 0.9) * rand) * HVL_params.y0;
        z1 = r1 + 1i*x1;
        z0 = r0 + 1i*x0;
        
        dist_km = ...
            FAULT_LOCATION_EXPRESSIONS_CORRECTED_NEW_2021...
            (Ubeg_PMU_set, Ibeg_PMU_set, Uend_PMU_set, Iend_PMU_set, ...
            z1, z0, y1, y0, len, sc_position, ...
            1:calc_settings.alg_count, sc_settings.faulted_phase); 

        delta_percent(idx_variation, :) = 100 * (dist_km - sc_position) / ...
            HVL_params.length; 
    end
    
    results_file_name = "delta_Monte_Carlo_line_params_variation.mat";
    results_path = fullfile(directories_settings.results_data, results_file_name);
    save(results_path, 'delta_percent');
end