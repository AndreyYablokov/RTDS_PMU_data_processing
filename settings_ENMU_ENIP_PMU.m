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
calc_settings.need_factors_ranges = false; % Необходим ли расчет диапазонов для факторов
calc_settings.threshold_val = 0.01; % Значение тока, по которому определяется начало выдачи тока RTDS
calc_settings.vectors_between_sc  = 15 * 100; % Интервал в кадрах между КЗ

% Варьируемые факторы
% calc_settings.factors = ["load" "sc_position" "sc_phase" "ef_mutual_angle" ...
%     "trans_resistance" "positive_seq" "zero_seq"];
calc_settings.factors = ["sc_position"];

% Нагрузка (в процентах от номинальной)
calc_settings.load.start_value = 0;
calc_settings.load.finish_value = 1;
calc_settings.load.step_change = 0.2;
calc_settings.load.range = [];
calc_settings.load.exp_num = 0;

% Местоположение КЗ (в долях от длины линии)
calc_settings.sc_position.start_value = 0;
calc_settings.sc_position.finish_value = 1;
calc_settings.sc_position.step_change = 0.01;
calc_settings.sc_position.range = [0.004545455 0.090909091 0.181818182 0.272727273 0.363636364 0.454545455 0.545454545 0.636363636 0.727272727 0.818181818 0.909090909 0.995454545];
calc_settings.sc_position.exp_num = 12;

% Фаза КЗ (в градусах)
calc_settings.sc_phase.start_value = 0;
calc_settings.sc_phase.finish_value = 350;
calc_settings.sc_phase.step_change = 10;
calc_settings.sc_phase.range = [];
calc_settings.sc_phase.exp_num = 0;

% Взаимный угол ЭДС (в градусах)
calc_settings.ef_mutual_angle.start_value = -20;
calc_settings.ef_mutual_angle.finish_value = 70;
calc_settings.ef_mutual_angle.step_change = 10;
calc_settings.ef_mutual_angle.range = [];
calc_settings.ef_mutual_angle.exp_num = 0;

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
calc_settings.trans_resistance.range = [];
calc_settings.trans_resistance.exp_num = 0;

% Отношение сопротивления ПП системы справа к системе слева
calc_settings.positive_seq.start_value = 0.1;
calc_settings.positive_seq.finish_value = 2;
calc_settings.positive_seq.step_change = 0.1;
calc_settings.positive_seq.range = [];
calc_settings.positive_seq.exp_num = 0;

% Отношение сопротивления НП системы справа к системе слева
calc_settings.zero_seq.start_value = 0.1;
calc_settings.zero_seq.finish_value = 4.0;
calc_settings.zero_seq.step_change = 0.1;
calc_settings.zero_seq.range = [];
calc_settings.zero_seq.exp_num = 0;

% Параметры КЗ
sc_settings.fault_duration = 5 * milliseconds(20); % Длительность КЗ
sc_settings.prefault_duration = 10 * milliseconds(20); % Длительность нормального режима перед КЗ
sc_settings.position_init = 110; % Местоположение КЗ
sc_settings.faulted_phase = 1; % Поврежденная фаза (1 - А, 2 - В,3 - С)

% Моменты возникновения КЗ для каждого их экспериментов
sc_settings.fault_inception_moment.load = [1 1 1 1 1 1];
sc_settings.fault_inception_moment.sc_position = 1;
sc_settings.fault_inception_moment.sc_phase = 1;
sc_settings.fault_inception_moment.ef_mutual_angle = [1 1 1 1 1 1 1 1 1 1];
sc_settings.fault_inception_moment.trans_resistance = 1;
sc_settings.fault_inception_moment.positive_seq = 1;
sc_settings.fault_inception_moment.zero_seq = 1;

% Настройки данных PMU
PMU_settings.need_change_direction_Iend = true;
PMU_settings.angle_in_degrees = true;
PMU_settings.column_count = 24;
PMU_settings.period = 2;
PMU_settings.Ibeg.phase_A.amp = 1;
PMU_settings.Ibeg.phase_A.angle = 2;
PMU_settings.Ibeg.phase_B.amp = 3;
PMU_settings.Ibeg.phase_B.angle = 4;
PMU_settings.Ibeg.phase_C.amp = 5;
PMU_settings.Ibeg.phase_C.angle = 6;
PMU_settings.Iend.phase_A.amp = 7;
PMU_settings.Iend.phase_A.angle = 8;
PMU_settings.Iend.phase_B.amp = 9;
PMU_settings.Iend.phase_B.angle = 10;
PMU_settings.Iend.phase_C.amp = 11;
PMU_settings.Iend.phase_C.angle = 12;
PMU_settings.Ubeg.phase_A.amp = 13;
PMU_settings.Ubeg.phase_A.angle = 14;
PMU_settings.Ubeg.phase_B.amp = 15;
PMU_settings.Ubeg.phase_B.angle = 16;
PMU_settings.Ubeg.phase_C.amp = 17;
PMU_settings.Ubeg.phase_C.angle = 18;
PMU_settings.Uend.phase_A.amp = 19;
PMU_settings.Uend.phase_A.angle = 20;
PMU_settings.Uend.phase_B.amp = 21;
PMU_settings.Uend.phase_B.angle = 22;
PMU_settings.Uend.phase_C.amp = 23;
PMU_settings.Uend.phase_C.angle = 24;
PMU_settings.names.ENIP = {
        'Status'...
        'VphsAan', 'VphsAmag', 'VphsBan', 'VphsBmag', 'VphsCan', 'VphsCmag'... 
        'IphsAan', 'IphsAmag', 'IphsBan', 'IphsBmag', 'IphsCan', 'IphsCmag'...
        'V0an', 'V0mag', 'V1an', 'V1mag', 'V2an', 'V2mag'...
        'I0an', 'I0mag', 'I1an', 'I1mag', 'I2an', 'I2mag'...
        'Frequency', 'dfdt'...
        'VrmsA', 'IrmsA','VrmsB', 'IrmsB','VrmsC', 'IrmsC'...
        'fA', 'fB', 'fC', 'df', 'EOF'
        };
PMU_settings.names.ENMU = {
        'Status'...
        'VphsAan', 'VphsAmag', 'VphsBan', 'VphsBmag', 'VphsCan', 'VphsCmag'... 
        'IphsAanCC', 'IphsAmagCC', 'IphsBanCC', 'IphsBmagCC', 'IphsCanCC', 'IphsCmagCC'...
        'IphsAan', 'IphsAmag', 'IphsBan', 'IphsBmag', 'IphsCan', 'IphsCmag'...
        'Frequency', 'dfdt', 'EOF'
        };
PMU_settings.PMU_data_beg_line = 'ENMU'; % PMU в начале линии
PMU_settings.PMU_data_end_line = 'ENIP'; % PMU в конце линии
PMU_settings.is_not_runtime_PMU = true; % Если данные не из runtime RTDS, то должно быть true
PMU_settings.PMU_shift_beg_line = 3.5 * milliseconds(20); % Сдвиг фильтра относительно мгновенных значений
PMU_settings.PMU_shift_end_line = 3.5 * milliseconds(20); % Сдвиг фильтра относительно мгновенных значений

% Настройки SV
SV_settings.delta_time = hours(3) + seconds(37);
SV_settings.names.SV_beg_line = {
    'IphsA', 'IphsB', 'IphsC', 'I0'...
    'VphsA', 'VphsB', 'VphsC', 'V0'...
    };
SV_settings.names.SV_end_line = {
    'IphsA', 'IphsB', 'IphsC', 'I0'...
    'VphsA', 'VphsB', 'VphsC', 'V0'...
    };

% Настройки сохранения графиков
graphs_settings.need_save = 1; % Логический флаг сохранения графиков (1 - да, 0 - нет)
graphs_settings.xlim_min = 1; % Минимум по оси X
graphs_settings.xlim_max = 20; % Максимум по оси X
graphs_settings.ylim_min = 0; % Минимум по оси Y
graphs_settings.ylim_max = 90; % Максимум по оси Y
graphs_settings.algs_for_plot = [6 11 13 17 18]; % Алгоритмы для которых строим графики
graphs_settings.marker_size = 20; % Размер маркеров
graphs_settings.markers = ["Square" "Diamond" "o" "*" "+"]; % Маркеры для каждого алгоритма
graphs_settings.result_directory = "graphs"; % Директория для сохранения графиков 

% Расположение исходных данных и результатов
directories_settings.PMU_data_MAT = 'PMU_data_MAT'; % Директория с данными PMU в MAT-файлах
directories_settings.PMU_data_CSV = 'PMU_data_CSV'; % Директория с данными PMU в CSV-файлах
directories_settings.SV_data_CSV = 'SV_data_CSV'; % Директория с SV в CSV-файлах
directories_settings.results_data = 'results'; % Директория с результатами расчета