package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

// Параметры задачи
const (
	A        = 0.0
	B        = 4.0
	M        = 11 // Количество точек таблицы
	N        = 3  // Количество неизвестных параметров (степень полинома + 1)
	NumCheck = 21 // Количество проверочных точек
)

// Точка для графика
type Point struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
}

// Данные для графиков
type ChartData struct {
	FuncPoints   []Point   `json:"func_points"`   // f(x)
	ApproxPoints []Point   `json:"approx_points"` // phi(x)
	ErrorPoints  []Point   `json:"error_points"`  // d(x)
	Coefficients []float64 `json:"coefficients"`  // c0, c1, c2
}

// Исходная функция f(x) = sin^2(x) - x/5
func f(x float64) float64 {
	return math.Pow(math.Sin(x), 2) - x/5.0
}

// Решение СЛАУ методом Гаусса с выбором главного элемента
func gaussSolve(matrix [][]float64, rhs []float64) []float64 {
	n := len(rhs)
	// Создаём расширенную матрицу
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n+1)
		copy(a[i], matrix[i])
		a[i][n] = rhs[i]
	}

	// Прямой ход
	for i := 0; i < n; i++ {
		// Поиск главного элемента
		maxRow := i
		for k := i + 1; k < n; k++ {
			if math.Abs(a[k][i]) > math.Abs(a[maxRow][i]) {
				maxRow = k
			}
		}
		a[i], a[maxRow] = a[maxRow], a[i]

		// Обнуление элементов ниже главной диагонали
		for k := i + 1; k < n; k++ {
			factor := a[k][i] / a[i][i]
			for j := i; j <= n; j++ {
				a[k][j] -= factor * a[i][j]
			}
		}
	}

	// Обратный ход
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := a[i][n]
		for j := i + 1; j < n; j++ {
			sum -= a[i][j] * x[j]
		}
		x[i] = sum / a[i][i]
	}
	return x
}

func main() {
	fmt.Println("==================================================")
	fmt.Println("Аппроксимация функции f(x) = sin^2(x) - x/5")
	fmt.Println("Метод наименьших квадратов (PNS), полином 2-й степени")
	fmt.Println("==================================================")

	// 1. Построение таблицы (xi, yi)
	fmt.Println("\n--- Таблица исходных данных ---")
	fmt.Printf("%-5s | %-10s | %-12s\n", "i", "x_i", "y_i")
	fmt.Println("--------------------------------")

	xs := make([]float64, M)
	ys := make([]float64, M)

	step := (B - A) / float64(M-1)
	for i := 0; i < M; i++ {
		x := A + float64(i)*step
		y := f(x)
		xs[i] = x
		ys[i] = y
		fmt.Printf("%-5d | %-10.4f | %-12.6f\n", i+1, x, y)
	}

	// 2. Формирование системы нормальных уравнений
	// Нам нужны суммы: sum(x^0), sum(x^1), ..., sum(x^4) и sum(y*x^0), ..., sum(y*x^2)
	sumsX := make([]float64, 2*N-1) // степени от 0 до 4
	sumsY := make([]float64, N)     // sum(y), sum(xy), sum(x^2y)

	for i := 0; i < M; i++ {
		x := xs[i]
		y := ys[i]

		// Вычисление степеней x
		pow := 1.0
		for k := 0; k < 2*N-1; k++ {
			sumsX[k] += pow
			pow *= x
		}

		// Вычисление сумм для правой части
		pow = 1.0
		for k := 0; k < N; k++ {
			sumsY[k] += y * pow
			pow *= x
		}
	}

	// Матрица системы (3x3)
	matrix := make([][]float64, N)
	for i := 0; i < N; i++ {
		matrix[i] = make([]float64, N)
		for j := 0; j < N; j++ {
			matrix[i][j] = sumsX[i+j]
		}
	}

	// Решение системы
	coeffs := gaussSolve(matrix, sumsY)

	fmt.Println("\n--- Коэффициенты аппроксимации ---")
	fmt.Printf("c0 = %.6f\n", coeffs[0])
	fmt.Printf("c1 = %.6f\n", coeffs[1])
	fmt.Printf("c2 = %.6f\n", coeffs[2])
	fmt.Printf("phi(x) = %.6f + %.6f*x + %.6f*x^2\n", coeffs[0], coeffs[1], coeffs[2])

	// 3. Вычисление в 21 точке и подготовка данных для графиков
	fmt.Println("\n--- Проверочные точки (21 точка) ---")
	fmt.Printf("%-5s | %-10s | %-12s | %-12s | %-12s\n", "j", "x_j", "f(x_j)", "phi(x_j)", "d(x_j)")
	fmt.Println("------------------------------------------------------------------")

	chartData := ChartData{
		Coefficients: coeffs,
	}

	stepCheck := (B - A) / float64(NumCheck-1)
	for j := 0; j < NumCheck; j++ {
		x := A + float64(j)*stepCheck
		fVal := f(x)
		phiVal := coeffs[0] + coeffs[1]*x + coeffs[2]*x*x
		errVal := fVal - phiVal

		fmt.Printf("%-5d | %-10.4f | %-12.6f | %-12.6f | %-12.6f\n", j+1, x, fVal, phiVal, errVal)

		chartData.FuncPoints = append(chartData.FuncPoints, Point{X: x, Y: fVal})
		chartData.ApproxPoints = append(chartData.ApproxPoints, Point{X: x, Y: phiVal})
		chartData.ErrorPoints = append(chartData.ErrorPoints, Point{X: x, Y: errVal})
	}

	// 4. Анализ погрешности
	maxErr := 0.0
	for _, p := range chartData.ErrorPoints {
		if math.Abs(p.Y) > maxErr {
			maxErr = math.Abs(p.Y)
		}
	}
	fmt.Println("------------------------------------------------------------------")
	fmt.Printf("Максимальная погрешность: %.6f\n", maxErr)

	// 5. Сохранение данных для графиков
	jsonData, err := json.MarshalIndent(chartData, "", "  ")
	if err != nil {
		fmt.Println("Ошибка сериализации JSON:", err)
		return
	}

	err = os.WriteFile("chart_data.json", jsonData, 0644)
	if err != nil {
		fmt.Println("Ошибка записи файла:", err)
	} else {
		fmt.Println("\nДанные для графиков сохранены в chart_data.json")
	}
}
