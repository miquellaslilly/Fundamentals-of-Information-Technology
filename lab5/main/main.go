package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

const (
	A       = -9.0
	B       = -1.0
	M       = 80   // количество разбиений таблицы
	Eps     = 1e-4 // базовая точность
	Q       = 4.0  // коэффициент дробления шага
	MaxIter = 100000
)

// Счётчик вызовов функции
var funcCalls int

// Исходная функция f(x) = 20x*sin^2(x) - x^2
func f(x float64) float64 {
	funcCalls++
	s := math.Sin(x)
	return 20*x*s*s - x*x
}

// Итерация MPP (для истории сходимости)
type Iteration struct {
	K  int     `json:"k"`
	X  float64 `json:"x"`
	H  float64 `json:"h"`
	Fx float64 `json:"fx"`
}

// Результат поиска одного минимума
type MinResult struct {
	X0         float64     `json:"x0"`
	H0         float64     `json:"h0"`
	Eps        float64     `json:"eps"`
	Minimum    float64     `json:"minimum"`
	Fmin       float64     `json:"fmin"`
	Iterations int         `json:"iterations"`
	Calls      int         `json:"calls"` // число вызовов f(x)
	History    []Iteration `json:"history"`
}

// Метод MPP (поразрядный поиск)
func mpp(x0, h0, eps float64) MinResult {
	res := MinResult{X0: x0, H0: h0, Eps: eps}
	x := x0
	h := h0
	k := 0
	funcCalls = 0

	res.History = append(res.History, Iteration{K: 0, X: x, H: h, Fx: f(x)})

	for h > eps && k < MaxIter {
		k++
		// Пробуем шаг вправо
		f1 := f(x)
		f2 := f(x + h)

		if f2 < f1 {
			// функция убывает — двигаемся вправо
			x = x + h
		} else {
			// функция возрастает — уменьшаем шаг
			h = h / Q
		}
		res.History = append(res.History, Iteration{K: k, X: x, H: h, Fx: f(x)})
	}

	res.Minimum = x
	res.Fmin = f(x)
	res.Iterations = k
	res.Calls = funcCalls
	return res
}

// Вывод таблицы значений функции
func printTable(a, b float64, m int) {
	h := (b - a) / float64(m)
	fmt.Printf("\n--- Таблица значений f(x) на [%.1f, %.1f], шаг h = %.4f ---\n", a, b, h)
	fmt.Printf("%-6s | %-10s | %-12s\n", "i", "x", "f(x)")
	fmt.Println("----------------------------------")
	for i := 0; i <= m; i++ {
		x := a + float64(i)*h
		marker := ""
		if i > 0 {
			xPrev := a + float64(i-1)*h
			if f(x) < f(xPrev) {
				// локальный рост? проверим следующий
			}
			_ = xPrev
		}
		// Отмечаем смену знака производной — грубое приближение минимума
		if i > 0 && i < m {
			xPrev := a + float64(i-1)*h
			xNext := a + float64(i+1)*h
			if f(x) < f(xPrev) && f(x) < f(xNext) {
				marker = "  <-- локальный минимум"
			}
		}
		fmt.Printf("%-6d | %-10.4f | %-12.6f%s\n", i, x, f(x), marker)
	}
}

// Автопоиск приблизительных локальных минимумов (для x0)
func findLocalMinima(a, b float64, m int) []float64 {
	h := (b - a) / float64(m)
	var mins []float64
	for i := 1; i < m; i++ {
		xPrev := a + float64(i-1)*h
		x := a + float64(i)*h
		xNext := a + float64(i+1)*h
		if f(x) < f(xPrev) && f(x) < f(xNext) {
			mins = append(mins, x)
		}
	}
	return mins
}

func main() {
	fmt.Println("=====================================================================")
	fmt.Println("Поиск локальных минимумов функции f(x) = 20x*sin^2(x) - x^2")
	fmt.Println("Метод MPP (поразрядный поиск), интервал [-9, -1]")
	fmt.Println("=====================================================================")

	// 1. Таблица значений
	printTable(A, B, M)

	// 2. Грубый поиск локальных минимумов (для начальных x0)
	seedMins := findLocalMinima(A, B, M)
	fmt.Printf("\nНайдено (грубо) локальных минимумов: %d\n", len(seedMins))
	for i, x := range seedMins {
		fmt.Printf("  %d: x0 ≈ %.4f, f = %.6f\n", i+1, x, f(x))
	}

	if len(seedMins) == 0 {
		fmt.Println("Локальные минимумы не найдены грубым поиском.")
		return
	}

	// 3. Уточнение каждого минимума методом MPP при разных eps
	epsValues := []float64{1e-2, 1e-3, 1e-4, 1e-5}
	h0 := (B - A) / 4.0 // начальный шаг

	fmt.Println("\n=====================================================================")
	fmt.Println("Уточнение минимумов методом MPP для разных eps")
	fmt.Println("=====================================================================")

	type EpsSummary struct {
		Eps        float64     `json:"eps"`
		Results    []MinResult `json:"results"`
		TotalCalls int         `json:"total_calls"`
	}
	var allSummaries []EpsSummary

	for _, eps := range epsValues {
		fmt.Printf("\n>>> Точность eps = %.0e\n", eps)
		var summary EpsSummary
		summary.Eps = eps
		totalCalls := 0

		for i, x0 := range seedMins {
			res := mpp(x0, h0, eps)
			summary.Results = append(summary.Results, res)
			totalCalls += res.Calls

			fmt.Printf("  Минимум %d: x* = %.8f, f(x*) = %.8f, итераций = %d, вызовов f = %d\n",
				i+1, res.Minimum, res.Fmin, res.Iterations, res.Calls)
		}
		summary.TotalCalls = totalCalls
		allSummaries = append(allSummaries, summary)
	}

	// 4. Сохранение данных для графиков
	outData := struct {
		Summaries []EpsSummary `json:"summaries"`
	}{Summaries: allSummaries}

	jsonData, err := json.MarshalIndent(outData, "", "  ")
	if err != nil {
		fmt.Println("Ошибка сериализации:", err)
		return
	}
	err = os.WriteFile("minima_data.json", jsonData, 0644)
	if err != nil {
		fmt.Println("Ошибка записи:", err)
	} else {
		fmt.Println("\nДанные сохранены в minima_data.json")
	}
}
