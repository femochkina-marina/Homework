# Оптимизация траектории полёта Ту-134

Программа рассчитывает оптимальную траекторию полёта самолёта.

## Что здесь есть

- `main.cpp` - расчёт траектории
- `plot_trajectory.cpp` - построение графиков

## Как запустить онлайн

### OnlineGDB
1. Зайти на https://www.onlinegdb.com/online_c++_compiler
2. Скопировать код из `main.cpp`
3. Нажать Run

### Replit
1. Зайти на https://replit.com
2. Создать C++ проект
3. Вставить код
4. Нажать Run

## Как установить на компьютер (Windows)

### Установка компилятора

**MinGW**
1. Скачать с https://sourceforge.net/projects/mingw/
2. Установить с галочкой mingw32-gcc-g++
3. Добавить в PATH папку C:\MinGW\bin

**Visual Studio**
1. Скачать с https://visualstudio.microsoft.com/ru/downloads/
2. Выбрать Desktop development with C++

### Компиляция и запуск

Открыть командную строку в папке с файлами:

```
g++ main.cpp -o main.exe
main.exe
```

Для графиков нужен gnuplot:

**Установка gnuplot:**
- Windows: скачать с http://www.gnuplot.info/download.html

```
g++ plot_trajectory.cpp -o plot.exe
plot.exe
```

## Что получится

После запуска main.exe создаются файлы:
- trajectory_vh.csv
- edges_A.csv
- edges_B.csv  
- edges_C.csv

После запуска plot.exe создаётся:
- trajectory_plots.png

## Режимы работы

```
main.exe vh fuel
main.exe vh time
main.exe vh mix 0.5
```

- fuel - экономия топлива
- time - быстрый полёт
- mix - комбинация
