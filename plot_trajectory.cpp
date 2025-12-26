#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>

using namespace std;

struct TrajectoryPoint {
    int i;
    char move;
    double V, H, t, m, alpha_deg, throttle, Theta_deg, dt, P;
};

vector<TrajectoryPoint> readCSV(const string& filename) {
    vector<TrajectoryPoint> data;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: файл " << filename << " не найден!" << endl;
        cerr << "Сначала запустите программу main для генерации данных." << endl;
        return data;
    }

    string line;
    getline(file, line);
    
    while (getline(file, line)) {
        TrajectoryPoint p;
        stringstream ss(line);
        string token;
        
        getline(ss, token, ','); p.i = stoi(token);
        getline(ss, token, ','); p.move = token.empty() ? '?' : token[0];
        getline(ss, token, ','); p.V = stod(token);
        getline(ss, token, ','); p.H = stod(token);
        getline(ss, token, ','); p.t = stod(token);
        getline(ss, token, ','); p.m = stod(token);
        getline(ss, token, ','); p.alpha_deg = stod(token);
        getline(ss, token, ','); p.throttle = stod(token);
        getline(ss, token, ','); p.Theta_deg = stod(token);
        getline(ss, token, ','); p.dt = stod(token);
        getline(ss, token, ','); p.P = stod(token);
        
        data.push_back(p);
    }
    
    file.close();
    return data;
}

void generateGnuplotScript(const vector<TrajectoryPoint>& data, const string& output_file) {
    if (data.empty()) return;
    
    ofstream gp("plot_script.gp");
    
    gp << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,10'\n";
    gp << "set output '" << output_file << "'\n";
    gp << "set encoding utf8\n";
    gp << "set multiplot layout 3,3 title 'Анализ траектории полета самолета Ту-134' font 'Arial,16'\n\n";
    
    ofstream tmpV("tmp_V.dat"), tmpH("tmp_H.dat"), tmpVH("tmp_VH.dat");
    ofstream tmpM("tmp_m.dat"), tmpF("tmp_fuel.dat"), tmpAlpha("tmp_alpha.dat");
    ofstream tmpThr("tmp_thr.dat"), tmpTheta("tmp_theta.dat"), tmpP("tmp_P.dat");
    
    double m0 = data[0].m;
    map<char, int> move_counts;
    
    for (const auto& p : data) {
        tmpVH << p.V << " " << p.H << " " << p.t << "\n";
        tmpV << p.t << " " << p.V << "\n";
        tmpH << p.t << " " << p.H << "\n";
        tmpM << p.t << " " << p.m << "\n";
        tmpF << p.t << " " << (m0 - p.m) << "\n";
        tmpAlpha << p.t << " " << p.alpha_deg << "\n";
        tmpThr << p.t << " " << (p.throttle * 100) << "\n";
        tmpTheta << p.t << " " << p.Theta_deg << "\n";
        tmpP << p.t << " " << (p.P / 1000.0) << "\n";
        if (p.move != '?') move_counts[p.move]++;
    }
    
    tmpV.close(); tmpH.close(); tmpVH.close(); tmpM.close(); tmpF.close();
    tmpAlpha.close(); tmpThr.close(); tmpTheta.close(); tmpP.close();
    
    gp << "set grid\n";
    gp << "set style line 1 lc rgb '#0060ad' lt 1 lw 2\n";
    gp << "set style line 2 lc rgb '#dd181f' lt 1 lw 2\n\n";
    
    gp << "set title 'Траектория в плоскости V-H' font 'Arial,12'\n";
    gp << "set xlabel 'Скорость V, м/с' font 'Arial,11'\n";
    gp << "set ylabel 'Высота H, м' font 'Arial,11'\n";
    gp << "plot 'tmp_VH.dat' u 1:2:3 with points pt 7 ps 0.5 lc palette notitle, \\\n";
    gp << "     'tmp_VH.dat' u 1:2 every ::0::0 with points pt 7 ps 2 lc rgb 'green' title 'Старт', \\\n";
    gp << "     'tmp_VH.dat' u 1:2 every ::" << (data.size()-1) << "::" << (data.size()-1) << " with points pt 5 ps 2 lc rgb 'red' title 'Финиш'\n\n";
    
    gp << "set title 'Скорость от времени' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Скорость V, м/с' font 'Arial,11'\n";
    gp << "plot 'tmp_V.dat' with lines ls 1 title 'Скорость'\n\n";
    
    gp << "set title 'Высота от времени' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Высота H, м' font 'Arial,11'\n";
    gp << "plot 'tmp_H.dat' with lines lc rgb '#00AA00' lw 2 title 'Высота'\n\n";
    
    gp << "set title 'Масса от времени' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Масса m, кг' font 'Arial,11'\n";
    gp << "set y2label 'Расход топлива, кг' font 'Arial,10' tc rgb 'orange'\n";
    gp << "set y2tics\n";
    gp << "plot 'tmp_m.dat' with lines ls 2 title 'Масса самолета', \\\n";
    gp << "     'tmp_fuel.dat' axes x1y2 with lines lc rgb 'orange' lw 2 dt 2 title 'Расход топлива'\n";
    gp << "unset y2label\n";
    gp << "unset y2tics\n\n";
    
    gp << "set title 'Угол атаки от времени' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Угол атаки α, град' font 'Arial,11'\n";
    gp << "plot 'tmp_alpha.dat' with lines lc rgb '#AA00AA' lw 2 title 'Угол атаки α', \\\n";
    gp << "     0 with lines lc rgb 'black' dt 2 notitle\n\n";
    
    gp << "set title 'Режим работы двигателя' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Режим двигателя, %' font 'Arial,11'\n";
    gp << "plot 'tmp_thr.dat' with lines lc rgb 'cyan' lw 2 title 'Режим двигателя', \\\n";
    gp << "     100 with lines lc rgb 'red' dt 2 title '100%'\n\n";
    
    gp << "set title 'Угол наклона траектории' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Угол траектории Θ, град' font 'Arial,11'\n";
    gp << "plot 'tmp_theta.dat' with lines lc rgb '#CCCC00' lw 2 title 'Угол Θ', \\\n";
    gp << "     0 with lines lc rgb 'black' dt 2 notitle\n\n";
    
    gp << "set title 'Тяга двигателя от времени' font 'Arial,12'\n";
    gp << "set xlabel 'Время t, с' font 'Arial,11'\n";
    gp << "set ylabel 'Тяга P, кН' font 'Arial,11'\n";
    gp << "plot 'tmp_P.dat' with lines lc rgb 'black' lw 2 title 'Тяга P'\n\n";
    
    gp << "set title 'Распределение типов маневров' font 'Arial,12'\n";
    gp << "set xlabel ' '\n";
    gp << "set ylabel 'Количество' font 'Arial,11'\n";
    gp << "set style data histograms\n";
    gp << "set style fill solid border -1\n";
    gp << "set xtics ('Тип A\\n(изм. скорости)' 0, 'Тип B\\n(изм. высоты)' 1, 'Тип C\\n(комбинир.)' 2) rotate by -15\n";
    gp << "plot '-' using 2:xtic(1) linecolor rgb variable notitle\n";
    gp << "A " << move_counts['A'] << " 0x87CEEB\n";
    gp << "B " << move_counts['B'] << " 0xF08080\n";
    gp << "C " << move_counts['C'] << " 0x90EE90\n";
    gp << "e\n\n";
    
    gp << "unset multiplot\n";
    
    gp.close();
    
    double total_time = data.back().t;
    double fuel_used = m0 - data.back().m;
    double v_start = data[0].V, v_end = data.back().V;
    double h_start = data[0].H, h_end = data.back().H;
    
    cout << string(60, '=') << endl;
    cout << "ВИЗУАЛИЗАЦИЯ ТРАЕКТОРИИ ПОЛЕТА" << endl;
    cout << string(60, '=') << endl;
    cout << "\nРЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ ТРАЕКТОРИИ Ту-134" << endl;
    cout << string(60, '-') << endl;
    cout << "Начальные условия:  V₀ = " << v_start << " м/с,  H₀ = " << h_start << " м" << endl;
    cout << "Конечные условия:   V_f = " << v_end << " м/с,  H_f = " << h_end << " м" << endl;
    cout << string(60, '-') << endl;
    cout << "Время полета:       " << total_time << " с  (" << (total_time/60.0) << " мин)" << endl;
    cout << "Расход топлива:     " << fuel_used << " кг" << endl;
    cout << "Начальная масса:    " << m0 << " кг" << endl;
    cout << "Конечная масса:     " << data.back().m << " кг" << endl;
    cout << "Количество точек:   " << data.size() << endl;
    cout << string(60, '=') << endl;
}

int main(int argc, char** argv) {
    string csv_file = (argc > 1) ? argv[1] : "trajectory_vh.csv";
    string output_file = "trajectory_plots.png";
    
    cout << "Чтение данных из: " << csv_file << endl;
    
    auto data = readCSV(csv_file);
    if (data.empty()) return 1;
    
    generateGnuplotScript(data, output_file);
    
    cout << "\nЗапуск gnuplot для генерации графиков..." << endl;
    int result = system("gnuplot plot_script.gp");
    
    if (result == 0) {
        cout << "✓ Графики сохранены в файл: " << output_file << endl;
        
        system("rm -f tmp_*.dat plot_script.gp");
        
        cout << "\nГотово!" << endl;
    } else {
        cerr << "Ошибка: не удалось запустить gnuplot." << endl;
        cerr << "Убедитесь, что gnuplot установлен: brew install gnuplot (macOS)" << endl;
        return 1;
    }
    
    return 0;
}

