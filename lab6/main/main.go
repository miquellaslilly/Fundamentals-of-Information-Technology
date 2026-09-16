package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

const (
	X0  = 1.0
	X1  = 4.0
	U10 = 2.0    // u1(1) = 2a, a=1
	U20 = math.E // u2(1) = e^1
	Eps = 1e-4
	H1  = (X1 - X0) / 10.0 // = 0.3
)

// Точное решение
func u1Exact(x float64) float64 { return 2 * x }
func u2Exact(x float64) float64 { return math.Exp(x) }

// Правая часть системы: F(x, u1, u2) = (f1, f2)
func FPR(x, u1, u2 float64) (float64, float64) {
	f1 := u1/x + u2 - math.Exp(x)
	f2 := 2*x/u1 + u2*u2/math.Exp(x) - 1.0
	return f1, f2
}

// Одна точка для таблицы
type TableRow struct {
	X   float64 `json:"x"`
	U1  float64 `json:"u1"`
	U2  float64 `json:"u2"`
	U1e float64 `json:"u1_exact"`
	U2e float64 `json:"u2_exact"`
	D1  float64 `json:"d1"`
	D2  float64 `json:"d2"`
}

// Решение задачи одним запуском RK4
type Solution struct {
	H        float64    `json:"h"`
	Steps    int        `json:"steps"`
	Rows     []TableRow `json:"rows"`
	MaxD1End float64    `json:"max_d1_end"` // |d1| в x=b
	MaxD2End float64    `json:"max_d2_end"` // |d2| в x=b
}

// Метод Рунге-Кутты 4-го порядка (M8)
func solveRK4(h float64) Solution {
	steps := int(math.Round((X1 - X0) / h))
	sol := Solution{H: h, Steps: steps}

	x := X0
	u1 := U10
	u2 := U20

	// Записываем начальную точку
	sol.Rows = append(sol.Rows, makeRow(x, u1, u2))

	for i := 0; i < steps; i++ {
		// k1
		k1_1, k1_2 := FPR(x, u1, u2)
		// k2
		k2_1, k2_2 := FPR(x+h/2, u1+h/2*k1_1, u2+h/2*k1_2)
		// k3
		k3_1, k3_2 := FPR(x+h/2, u1+h/2*k2_1, u2+h/2*k2_2)
		// k4
		k4_1, k4_2 := FPR(x+h, u1+h*k3_1, u2+h*k3_2)

		u1 += h / 6.0 * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)
		u2 += h / 6.0 * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)
		x += h

		sol.Rows = append(sol.Rows, makeRow(x, u1, u2))
	}

	// Погрешность в конечной точке
	last := sol.Rows[len(sol.Rows)-1]
	sol.MaxD1End = math.Abs(last.D1)
	sol.MaxD2End = math.Abs(last.D2)
	return sol
}

func makeRow(x, u1, u2 float64) TableRow {
	u1e := u1Exact(x)
	u2e := u2Exact(x)
	return TableRow{
		X: x, U1: u1, U2: u2,
		U1e: u1e, U2e: u2e,
		D1: u1 - u1e, D2: u2 - u2e,
	}
}

// Вывод таблицы (аналог процедуры OUT)
func OUT(row TableRow) {
	fmt.Printf("x=%6.4f  u1=%10.6f  u1*=%10.6f  d1=%10.2e   u2=%10.6f  u2*=%10.6f  d2=%10.2e\n",
		row.X, row.U1, row.U1e, row.D1, row.U2, row.U2e, row.D2)
}

func main() {
	fmt.Println("=====================================================================")
	fmt.Println("Решение задачи Коши для системы ОДУ методом Рунге-Кутты 4-го порядка")
	fmt.Println("u1' = u1/x + u2 - e^x")
	fmt.Println("u2' = 2x/u1 + u2^2/e^x - 1")
	fmt.Println("x in [1, 4], u1(1) = 2, u2(1) = e")
	fmt.Println("Точное решение: u1 = 2x, u2 = e^x")
	fmt.Println("=====================================================================")

	// Последовательность сгущающихся сеток
	var solutions []Solution
	h := H1
	for {
		sol := solveRK4(h)
		solutions = append(solutions, sol)

		fmt.Printf("\n--- h = %.6f  (шагов: %d) ---\n", h, sol.Steps)
		fmt.Println("  Таблица (последние 3 точки):")
		start := len(sol.Rows) - 3
		if start < 0 {
			start = 0
		}
		for _, row := range sol.Rows[start:] {
			OUT(row)
		}

		fmt.Printf("  Погрешность в x=4:  |d1| = %.3e, |d2| = %.3e\n",
			sol.MaxD1End, sol.MaxD2End)

		// Проверка условия остановки
		if sol.MaxD1End < Eps && sol.MaxD2End < Eps {
			fmt.Printf("\n>> Достигнута точность eps = %.0e при h = %.6f\n", Eps, h)
			break
		}

		h = h / 2.0
		if h < 1e-6 {
			fmt.Println(">> Достигнута предельная точность, дальнейшее уменьшение h нецелесообразно")
			break
		}
	}

	// Полная таблица для h = H1 (для графика)
	fmt.Println("\n=====================================================================")
	fmt.Println("Полная таблица решения при h = h1 = 0.3:")
	fmt.Println("=====================================================================")
	fmt.Printf("%-8s | %-10s | %-10s | %-10s | %-10s | %-10s\n",
		"x", "u1", "u1*", "d1", "u2", "u2*")
	fmt.Println("--------------------------------------------------------------------------")
	for _, row := range solutions[0].Rows {
		fmt.Printf("%-8.4f | %-10.6f | %-10.6f | %-10.2e | %-10.6f | %-10.6f\n",
			row.X, row.U1, row.U1e, row.D1, row.U2, row.U2e)
	}

	// Сохранение данных для графиков
	jsonData, err := json.MarshalIndent(solutions, "", "  ")
	if err != nil {
		fmt.Println("Ошибка сериализации:", err)
		return
	}
	err = os.WriteFile("ode_data.json", jsonData, 0644)
	if err != nil {
		fmt.Println("Ошибка записи:", err)
	} else {
		fmt.Println("\nДанные сохранены в ode_data.json")
	}
}
