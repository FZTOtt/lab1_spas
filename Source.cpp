#pragma warning(suppress : 4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>
#define OK 0
#define ERR -1


const double eps = 1e-15;  // машинный ноль
int width = 7;  // ширина пластины
int height = 3;  // высота пластины

// координаты центра окружности скругления

double D = 1;  // коэффициент теплопроводности
double delta = 0.25;  // шаг сетки по X и Y
double deltaT = 1;  // шаг моделирования в сек
int time1 = 30;  // время моделирования в сек
double q;  // D * deltaT / delta^2

double getTime() {
    return clock() / (CLOCKS_PER_SEC / 1000);
}

double* gauss(double** a, double* b, int n) {
    double* x = (double*)calloc(n, sizeof(double));
    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Находим максимальный элемент в столбце k
        double maxEl = a[k][k];
        int indMax = k;
        for (int i = k + 1; i < n; ++i) {
            if (fabs(maxEl) < fabs(a[i][k])) {
                maxEl = a[i][k];
                indMax = i;
            }
        }
        // Меняем местами
        if (k != indMax) {
            // Меняем строку с максимальным элементом местами с первой строкой
            double* tmpLine = a[k];
            a[k] = a[indMax];
            a[indMax] = tmpLine;
            // Меняем соотсветствующие элементы в векторе b
            double tmpEl = b[k];
            b[k] = b[indMax];
            b[indMax] = tmpEl;
        }
        if (fabs(a[k][k]) < eps) {
            printf("Решение получить невозможно из-за нулевого столбца ");
            printf("%d матрицы A\n", k);
            exit(-2);
        }
        // Вычитаем из строки i строку k умноженную на коэффициент coeff
        omp_set_num_threads(8);
        #pragma omp parallel for 
        for (int i = k + 1; i < n; ++i) {
            double coeff = a[i][k] / a[k][k];
            
            #pragma novector
            for (int j = k; j < n; ++j) {
                a[i][j] -= a[k][j] * coeff;
            }
            b[i] -= b[k] * coeff;
        }
    }
    // Делим каждую строку матрицы на соответствующий диагональный элемент
    for (int i = 0; i < n; ++i) {
        double coeff = a[i][i];
        for (int j = i; j < n; ++j) {
            if (fabs(a[i][j]) < eps) continue;
            a[i][j] /= coeff;
        }
        b[i] /= coeff;
    }
    // Обратный ход
    for (int i = n - 1; i >= 0; --i) {
        x[i] = a[i][i] * b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= a[i][j] * x[j];
        }
    }
    return x;
}


bool inPlace(int i, int j) {
    if (j == width - 1) {  // правая граница, которая тепловыделяющая
        return false;
    }
    return true;  // остальные границы и внутренние точки
}

bool inBorder(int i, int j) {
    return false;  
}

double getLambda(int i, int j) {
    return 0; 
}

double getMu(int i, int j) {
    return 0;  
}

void printTempMatrix(double** m, int height, int width) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {

            printf("%.1f\t", m[i][j]);

        }
        printf("\n");
    }
    printf("\n");
}

void addMatrixToFile(double** tempMatrix, FILE* file) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            fprintf(file, "%f %f %f\n", j * delta, (height - 1 - i) * delta, tempMatrix[i][j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n\n");

}


void gnuplot() {
    printf("gnuplotina");
    FILE* plot = fopen("plot.gpi", "w");
    char s[30];
    fprintf(plot, "set pm3d map\n");
    fprintf(plot, "set cbrange [0:250]\n");
    fprintf(plot, "set pm3d flush begin ftriangles scansforward interpolate 2,2\n");
    int K = (int)(time1 / deltaT);
    printf("%d\n", K);
    fprintf(plot, "do for [i=0:%d] {\n", K - 1);
    fprintf(plot, "\tsplot 'result.txt' index i using 1:2:3 with pm3d title '10 var'\n");
    sprintf(s, "\tpause 0.5\n");
    fprintf(plot, "%s", s);
    fprintf(plot, "}\n");
    fprintf(plot, "pause mouse key\n");
    fflush(plot);
    fclose(plot);
    system("gnuplot plot.gpi");
}

double** initTempMatrix() {
    // Выделение памяти под матрицу
    double** matrix = (double**)calloc(height, sizeof(double*));
    for (int i = 0; i < height; ++i) {
        matrix[i] = (double*)calloc(width, sizeof(double));
    }
    // Начальные значения температур
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if (j == width - 1) {
                matrix[i][j] = 0;  // правая граница (тепловыделяющая)
            }
            else {
                matrix[i][j] = 100;  // остальные границы и внутренние точки
            }
        }
    }

    return matrix;
}

void nextTmpMatrix(double** tempMatrix) {
    // Выделение памяти под матрицу
    double** matrix = (double**)calloc(width * height, sizeof(double*));
    for (int i = 0; i < width * height; ++i) {
        matrix[i] = (double*)calloc(width * height, sizeof(double));
    }
    // Выделение памяти под вектор b
    double* b = (double*)calloc(width * height, sizeof(double));


    // Заполнение матрицы
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            // Граничное условие 1го рода (сохранение температуры)
            if (i == 0 || j == 0 || i == height - 1 || j == width - 1) {
                matrix[i * width + j][i * width + j] = -1;
            }
            else {  // На границе скругления и за ее пределами
                if (!inPlace(i, j)) {
                    matrix[i * width + j][i * width + j] = -1;
                }
                else {  // Разностная схема уравнения теплопроводности
                    matrix[i * width + j][i * width + j] = -1;
                    // Если верхняя соседняя точка находится на скруглении
                    if (inBorder(i - 1, j)) {
                        double lambda = getLambda(i - 1, j);
                        matrix[i * width + j][(i - 1) * width + j] = 2 * q / (lambda
                            + 1);
                        matrix[i * width + j][i * width + j] -= 2 * q / lambda;
                        matrix[i * width + j][(i + 1) * width + j] = 2 * q / (lambda
                            * (lambda + 1));
                    }
                    else {
                        matrix[i * width + j][(i - 1) * width + j] = q;
                        matrix[i * width + j][i * width + j] -= 2 * q;
                        matrix[i * width + j][(i + 1) * width + j] = q;
                    }
                    // Если правая соседняя точка находится на скруглении
                    if (inBorder(i, j + 1)) {
                        double mu = getMu(i, j + 1);
                        matrix[i * width + j][i * width + j - 1] = 2 * q / (mu + 1);
                        matrix[i * width + j][i * width + j] -= 2 * q / mu;
                        matrix[i * width + j][i * width + j + 1] = 2 * q / (mu * (mu
                            + 1));
                    }
                    else {


                        matrix[i * width + j][i * width + j - 1] = q;
                        matrix[i * width + j][i * width + j] -= 2 * q;
                        matrix[i * width + j][i * width + j + 1] = q;

                    }
                }
            }
        }
    }
    // Заполнение вектора b
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            b[i * width + j] = -tempMatrix[i][j];
        }
    }
    // Решение полученной СЛАУ
    double* result = gauss(matrix, b, width * height);
    // Заполнение матрицы температурного распределения
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            tempMatrix[i][j] = result[i * width + j];
        }
    }
    free(result);
    // Освобождение памяти матрицы расчета
    for (int j = 0; j < width * height; ++j) {
        free(matrix[j]);
    }
    free(matrix);
    // Освобождение памяти вектора b
    free(b);
}

int main(int argc, char** argv) {
    double t_start = getTime();
    // Приведение размеров с учетом значения delta
    width = width / delta + 1;
    height = height / delta + 1;

    q = D * deltaT / pow(delta, 2);

    // Открытие файлов для записи результата
    FILE* file;
    if ((file = fopen("result.txt", "w")) == NULL) {  // открытие файла для записирезультата
        printf("Ошибка! Не удалось открыть файл для записи результатов.\n");
        return ERR;
    }
    // Инициализация матрицы начальных температур
    double** tempMatrix = initTempMatrix();
    addMatrixToFile(tempMatrix, file);

    // Вывод матрицы температур начального момента времени в консоль
    printTempMatrix(tempMatrix, height, width);


    // Основной цикл расчета
    for (int t = 0; t < time1 / deltaT; ++t) {
        printf("%d\n", t);
        nextTmpMatrix(tempMatrix);
        addMatrixToFile(tempMatrix, file);
    }

    // Вывод матрицы температур конечного момента времени в консоль
    printTempMatrix(tempMatrix, height, width);

    // Освобождение памяти матрицы температурных значений
    for (int j = 0; j < height; ++j) {
        free(tempMatrix[j]);
    }
    free(tempMatrix);

    double t_end = getTime();
    double t = t_end - t_start;

    gnuplot();

    printf("%f", t);

    return OK;
}
