clear
%% Исходные данные
% Параметры ВЛ
HVL_params.length = 220; % Длина ВЛ, км
HVL_params.r0 = 0.207482; % Погонное активное сопротивление нулевой последовательности ВЛ, Ом/км
HVL_params.r1 = 0.0253524; % Погонное активное сопротивление прямой последовательности ВЛ, Ом/км
HVL_params.x0 = 0.733642; % Погонное реактивное сопротивление нулевой последовательности ВЛ, Ом/км
HVL_params.x1 = 0.302264; % Погонное реактивное сопротивление прямой последовательности ВЛ, Ом/км
HVL_params.y0 = 2.63318e-06; % Погонная проводимость нулевой последовательности ВЛ, См/км
HVL_params.y1 = 3.80648e-06; % Погонная проводимость прямой последовательности ВЛ, См/км

% Параметры расчета
calc_settings.need_load_PMU_from_csv = true; % Необходимо загружать данные из csv или есть данные в MAT-файле
calc_settings.alg_num = 13; % Для какого алгоритма строим временные зависимости?
calc_settings.alg_count = 18; % Количество исследуемых алгоритмов ОМП
calc_settings.calc_period_num = 3; % Через сколько периодов выбираем комплект СВИ для расчета?

% Варьируемые факторы
% calc_settings.factors = ["load" "sc_position" "sc_phase" "ef_mutual_angle" ...
%     "trans_resistance" "positive_seq" "zero_seq"];
calc_settings.factors = ["sc_position"];

% Нагрузка (в процентах от номинальной)
calc_settings.load.start_value = 0;
calc_settings.load.finish_value = 1;
calc_settings.load.step_change = 0.2;

% Местоположение КЗ (в процентах от длины линии)
calc_settings.sc_position.start_value = 0;
calc_settings.sc_position.finish_value = 1;
calc_settings.sc_position.step_change = 0.1;

% Фаза КЗ (в градусах)
calc_settings.sc_phase.start_value = 0;
calc_settings.sc_phase.finish_value = 350;
calc_settings.sc_phase.step_change = 10;

% Взаимный угол ЭДС (в градусах)
calc_settings.ef_mutual_angle.start_value = -20;
calc_settings.ef_mutual_angle.finish_value = 70;
calc_settings.ef_mutual_angle.step_change = 10;

% Переходное сопротивление (в Ом)
calc_settings.trans_resistance.start_value_int1 = 0;
calc_settings.trans_resistance.finish_value_int1 = 1;
calc_settings.trans_resistance.step_change_int1 = 0.1;
calc_settings.trans_resistance.start_value_int2 = 2;
calc_settings.trans_resistance.finish_value_int2 = 10;
calc_settings.trans_resistance.step_change_int2 = 1;
calc_settings.trans_resistance.start_value_int3 = 20;
calc_settings.trans_resistance.finish_value_int3 = 50;
calc_settings.trans_resistance.step_change_int3 = 10;

% Отношение сопротивления ПП системы справа к системе слева
calc_settings.positive_seq.start_value = 0.1;
calc_settings.positive_seq.finish_value = 2;
calc_settings.positive_seq.step_change = 0.1;

% Отношение сопротивления НП системы справа к системе слева
calc_settings.zero_seq.start_value = 0.1;
calc_settings.zero_seq.finish_value = 4.0;
calc_settings.zero_seq.step_change = 0.1;

% Параметры КЗ
sc_settings.position_init = 110; % Местоположение КЗ
sc_settings.faulted_phase = 1; % Поврежденная фаза (1 - А, 2 - В,3 - С)

% Моменты возникновения КЗ для каждого их экспериментов
sc_settings.fault_inception_moment.load = [2095 2131 2164 2177 2184 2188];
sc_settings.fault_inception_moment.sc_position = 2188;
sc_settings.fault_inception_moment.sc_phase = 2188;
sc_settings.fault_inception_moment.ef_mutual_angle = [2188 2187 2347 2363 2360 2353 2350 2345 2340 2334];
sc_settings.fault_inception_moment.trans_resistance = 2188;
sc_settings.fault_inception_moment.positive_seq = 2188;
sc_settings.fault_inception_moment.zero_seq = 2188;

% Настройки данных PMU
PMU_settings.column_count = 31;
PMU_settings.period = 400;
PMU_settings.Ibeg.phase_A.amp = 3;
PMU_settings.Ibeg.phase_A.angle = 2;
PMU_settings.Ibeg.phase_B.amp = 5;
PMU_settings.Ibeg.phase_B.angle = 4;
PMU_settings.Ibeg.phase_C.amp = 7;
PMU_settings.Ibeg.phase_C.angle = 6;
PMU_settings.Iend.phase_A.amp = 18;
PMU_settings.Iend.phase_A.angle = 17;
PMU_settings.Iend.phase_B.amp = 20;
PMU_settings.Iend.phase_B.angle = 19;
PMU_settings.Iend.phase_C.amp = 22;
PMU_settings.Iend.phase_C.angle = 21;
PMU_settings.Ubeg.phase_A.amp = 12;
PMU_settings.Ubeg.phase_A.angle = 11;
PMU_settings.Ubeg.phase_B.amp = 14;
PMU_settings.Ubeg.phase_B.angle = 13;
PMU_settings.Ubeg.phase_C.amp = 16;
PMU_settings.Ubeg.phase_C.angle = 15;
PMU_settings.Uend.phase_A.amp = 27;
PMU_settings.Uend.phase_A.angle = 26;
PMU_settings.Uend.phase_B.amp = 29;
PMU_settings.Uend.phase_B.angle = 28;
PMU_settings.Uend.phase_C.amp = 31;
PMU_settings.Uend.phase_C.angle = 30;

% Настройки сохранения графиков
graphs_settings.need_save = 1; % Логический флаг сохранения графиков (1 - да, 0 - нет)
graphs_settings.xlim_min = 1800; % Минимум по оси X
graphs_settings.xlim_max = 5000; % Максимум по оси X
graphs_settings.ylim_min = 0; % Минимум по оси Y
graphs_settings.ylim_max = 90; % Максимум по оси Y
graphs_settings.algs_for_plot = [6 11 13 17 18]; % Алгоритмы для которых строим графики
graphs_settings.marker_size = 20; % Размер маркеров
graphs_settings.markers = ["Square" "Diamond" "o" "*" "+"]; % Маркеры для каждого алгоритма
graphs_settings.result_directory = "graphs"; % Директория для сохранения графиков 

% Расположение исходных данных и результатов
directories_settings.PMU_data_MAT = 'PMU_data_MAT'; % Директория с исходными данными в MAT-файлах
directories_settings.PMU_data_CSV = 'PMU_data_CSV'; % Директория с исходными данными в CSV-файлах
directories_settings.results_data = 'results'; % Директория с результатами расчета

%% Расчет недостающих данных
% Параметры ВЛ
HVL_params.z0 = HVL_params.r0 + HVL_params.x0 * 1i; % Погонное полное сопротивление нулевой последовательности ВЛ
HVL_params.z1 = HVL_params.r1 + HVL_params.x1 * 1i; % Погонное полное сопротивление прямой последовательности ВЛ

% Параметры расчета
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

%% Основная прграмма
for factor = calc_settings.factors

    %% Загрузка данных УСВИ
    if calc_settings.need_load_PMU_from_csv == true
        PMU_data_matrix = load_PMU_data_from_csv(...
            directories_settings, factor, calc_settings);
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
        delta_percent = Monte_Carlo_data_processing(...
            Ibeg, Iend, Ubeg, Uend, factor, PMU_settings.period, HVL_params, ...
            calc_settings, sc_settings, directories_settings);
        
        create_errors_graphs_Monte_Carlo(...
            delta_percent, graphs_settings)
    end
    
    %% Освобождение памяти от лишних переменных
    clearvars -except HVL_params calc_settings sc_settings PMU_settings ...
        graphs_settings graphs_settings directories_settings;
    
end

clear;

%% Функции

% Функция загрузки данных из CSV
function PMU_data_matrix = load_PMU_data_from_csv(directories_settings, factor, calc_settings)
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

% Функция выделения токов и напряжений из матрицы данных PMU
function [Ibeg, Iend, Ubeg, Uend] = get_IU_from_PMU(PMU, factor, PMU_settings, calc_settings)

    exp_count = getfield(calc_settings,factor).exp_num;   

    for idx_exp = 1:exp_count
        Ibeg.phase_A(:,idx_exp) = PMU(:,PMU_settings.Ibeg.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Ibeg.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Ibeg.phase_B(:,idx_exp) = PMU(:,PMU_settings.Ibeg.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Ibeg.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Ibeg.phase_C(:,idx_exp) = PMU(:,PMU_settings.Ibeg.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Ibeg.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));

        Iend.phase_A(:,idx_exp) = PMU(:,PMU_settings.Iend.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Iend.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Iend.phase_B(:,idx_exp) = PMU(:,PMU_settings.Iend.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Iend.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Iend.phase_C(:,idx_exp) = PMU(:,PMU_settings.Iend.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Iend.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
      
        Ubeg.phase_A(:,idx_exp) = PMU(:,PMU_settings.Ubeg.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Ubeg.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Ubeg.phase_B(:,idx_exp) = PMU(:,PMU_settings.Ubeg.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Ubeg.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Ubeg.phase_C(:,idx_exp) = PMU(:,PMU_settings.Ubeg.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Ubeg.phase_C.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));

        Uend.phase_A(:,idx_exp) = PMU(:,PMU_settings.Uend.phase_A.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Uend.phase_A.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Uend.phase_B(:,idx_exp) = PMU(:,PMU_settings.Uend.phase_B.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Uend.phase_B.angle + ...
            ((idx_exp-1)*PMU_settings.column_count))));
        
        Uend.phase_C(:,idx_exp) = PMU(:,PMU_settings.Uend.phase_C.amp + ...
            ((idx_exp-1)*PMU_settings.column_count)).*...
            exp(1j*(PMU(:,PMU_settings.Uend.phase_C.angle + ...
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

    sc_position = sc_settings.position_init;
    exp_count = getfield(calc_settings,factor).exp_num;
    
    delta_percent = zeros(exp_count, 20000);

    for idx_exp = 1:exp_count 
        for PMU_set = 1:20000         
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
            end
                
            dist_km = ...
                FAULT_LOCATION_EXPRESSIONS_CORRECTED_NEW_2021 ...
                (Ubeg_PMU_set, Ibeg_PMU_set, Uend_PMU_set, Iend_PMU_set, ...
                HVL_params.z1, HVL_params.z0, HVL_params.y1, HVL_params.y0, ...
                HVL_params.length, sc_position, calc_settings.alg_num, sc_settings.faulted_phase);

                delta_percent(idx_exp, PMU_set) = ...
                100 * ((dist_km - sc_position) ...
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
                saveas(gcf, results_path, 'emf' );
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
            saveas(gcf, results_path, 'emf' );
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

% Функция расчета погрешностей алгоритмов по методу Монте-Карло
function delta_percent = Monte_Carlo_data_processing(...
    Ibeg, Iend, Ubeg, Uend, factor, period, HVL_params, ...
    calc_settings, sc_settings, directories_settings)
    
    sc_position = sc_settings.position_init;
    fault_inception_moment = getfield(sc_settings.fault_inception_moment,factor);
    PMU_set = fault_inception_moment + period;
    
    if factor ~= "load" && factor ~= "sc_phase" && factor ~= "ef_mutual_angle"
        PMU_set = fault_inception_moment + 3*period;
    else
        PMU_set = fault_inception_moment(1) + 3*period;
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
    
    results_file_name = "delta_Monte_Carlo.mat";
    results_path = fullfile(directories_settings.results_data, results_file_name);
    save(results_path, 'delta_percent');
end

% Функция построения графиков погрешностей, определенных по методу
% Монте-Карло
function create_errors_graphs_Monte_Carlo(...
    delta_percent, graphs_settings)

    for idx_alg = graphs_settings.algs_for_plot
        fig1 = figure;
        probplot('normal', real(delta_percent(:, idx_alg)));
        title(strcat(num2str(idx_alg),' алг'));
        
        if graphs_settings.need_save == 1
            graph_file_name = strcat("Monte_Carlo_res_prob_alg", num2str(idx_alg));
            results_path = fullfile(graphs_settings.result_directory, graph_file_name);
            saveas(gcf, results_path, 'png' );
            saveas(gcf, results_path, 'fig' );
            saveas(gcf, results_path, 'emf' );
        end 
        
        close (fig1);
        
        fig2 = figure;
        histfit(real(delta_percent(:, idx_alg)));
        hold on;
        xline(mean(real(delta_percent(:, idx_alg))),'r');
        hold on;
        xline(mean(real(delta_percent(:, idx_alg))) - ...
            3*std(real(delta_percent(:, idx_alg))),'r')
        hold on;
        xline(mean(real(delta_percent(:, idx_alg))) + ...
            3*std(real(delta_percent(:, idx_alg))),'r')
        title(strcat(num2str(idx_alg),' алг'));
        
        if graphs_settings.need_save == 1
            graph_file_name = strcat("Monte_Carlo_res_hist_alg", num2str(idx_alg));
            results_path = fullfile(graphs_settings.result_directory, graph_file_name);
            saveas(gcf, results_path, 'png' );
            saveas(gcf, results_path, 'fig' );
            saveas(gcf, results_path, 'emf' );
        end 
        
        close (fig2);
    end
end
