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
calc_settings.alg_num = 17; % Для какого алгоритма строим временные зависимости?
calc_settings.alg_count = 18; % Количество исследуемых алгоритмов ОМП
calc_settings.calc_period_num = 3; % Через сколько периодов выбираем комплект СВИ для расчета?
calc_settings.need_factors_ranges = true; % Необходим ли расчет диапазонов для факторов

% Варьируемые факторы
% calc_settings.factors = ["load" "sc_position" "sc_phase" "ef_mutual_angle" ...
%     "trans_resistance" "positive_seq" "zero_seq" "Monte_Carlo"];
calc_settings.factors = ["sc_position"];

% Нагрузка (в процентах от номинальной)
calc_settings.load.start_value = 0;
calc_settings.load.finish_value = 1;
calc_settings.load.step_change = 0.2;
calc_settings.load.range = [];
calc_settings.load.exp_num = 0;

% Местоположение КЗ (в процентах от длины линии)
calc_settings.sc_position.start_value = 0;
calc_settings.sc_position.finish_value = 1;
calc_settings.sc_position.step_change = 0.01;
calc_settings.sc_position.range = [];
calc_settings.sc_position.exp_num = 0;

% Фаза КЗ (в градусах)
calc_settings.sc_phase.start_value = 0;
calc_settings.sc_phase.finish_value = 360;
calc_settings.sc_phase.step_change = 10;
calc_settings.sc_phase.range = [];
calc_settings.sc_phase.exp_num = 0;

% Взаимный угол ЭДС (в градусах)
calc_settings.ef_mutual_angle.start_value = -20;
calc_settings.ef_mutual_angle.finish_value = 90;
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

% Количество экспериментов метода Монте-Карло
calc_settings.Monte_Carlo.exp_num = 0;
calc_settings.Monte_Carlo.range = 1:1:calc_settings.Monte_Carlo.exp_num;

% Параметры КЗ
sc_settings.position_init = 110; % Местоположение КЗ
sc_settings.faulted_phase = 1; % Поврежденная фаза (1 - А, 2 - В,3 - С)

% Моменты возникновения КЗ для каждого их экспериментов
sc_settings.fault_inception_moment.load = [1000 1000 1000 1000 1000 1000];
sc_settings.fault_inception_moment.sc_position = 1000;
sc_settings.fault_inception_moment.sc_phase = 1000;
sc_settings.fault_inception_moment.ef_mutual_angle = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
sc_settings.fault_inception_moment.trans_resistance = 1000;
sc_settings.fault_inception_moment.positive_seq = 1000;
sc_settings.fault_inception_moment.zero_seq = 1000;
sc_settings.fault_inception_moment.Monte_Carlo = 1000;

% Настройки данных PMU
PMU_settings.need_change_direction_Iend = true;
PMU_settings.angle_in_degrees = false;
PMU_settings.RC_data = true;
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

PMU_settings.is_not_runtime_PMU = false; % Если данные не из runtime RTDS, то должно быть true

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