package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

const (
	A   = -4.0
	B   = 0.0
	M   = 40   // количество разбиений для таблицы
	Eps = 1e-4 // точность
)

// Исходная функция f(x) = sin^2(x) - x/5 - 1
func f(x float64) float64 {
	return math.Pow(math.Sin(x), 2) - x/5.0 - 1
}

// Производная f'(x) = sin(2x) - 0.2
func df(x float64) float64 {
	return math.Sin(2*x) - 0.2
}

// Итерация метода Ньютона (MN)
type Iteration struct {
	K   int     `json:"k"`
	Xk  float64 `json:"xk"`
	Fxk float64 `json:"fxk"`
	Dx  float64 `json:"dx"` // d_k = x_k - x_{k-1}
}

// Результат поиска одного корня
type RootResult struct {
	X0         float64     `json:"x0"`         // начальное приближение
	Root       float64     `json:"root"`       // найденный корень
	Iterations int         `json:"iterations"` // количество итераций
	History    []Iteration `json:"history"`    // история итераций
	Froot      float64     `json:"froot"`      // f(root)
}

// Метод Ньютона
func newton(x0 float64, eps float64, maxIter int) RootResult {
	res := RootResult{X0: x0}
	x := x0
	xPrev := x0
	k := 0

	res.History = append(res.History, Iteration{K: 0, Xk: x, Fxk: f(x), Dx: 0})

	for k = 1; k <= maxIter; k++ {
		fx := f(x)
		dfx := df(x)
		if math.Abs(dfx) < 1e-15 {
			fmt.Println("  Производная близка к нулю, метод Ньютона расходится")
			break
		}
		xNew := x - fx/dfx
		dx := xNew - x
		res.History = append(res.History, Iteration{K: k, Xk: xNew, Fxk: f(xNew), Dx: dx})

		xPrev = x
		x = xNew

		if math.Abs(dx) < eps || math.Abs(f(x)) < eps {
			break
		}
	}

	res.Root = x
	res.Iterations = k
	res.Froot = f(x)
	_ = xPrev
	return res
}

// Вывод таблицы значений
func printTable(a, b float64, m int) {
	h := (b - a) / float64(m)
	fmt.Printf("\n--- Таблица значений f(x) на [%.2f, %.2f], шаг h = %.4f ---\n", a, b, h)
	fmt.Printf("%-6s | %-10s | %-12s\n", "i", "x", "f(x)")
	fmt.Println("----------------------------------")
	for i := 0; i <= m; i++ {
		x := a + float64(i)*h
		fmt.Printf("%-6d | %-10.4f | %-12.6f\n", i, x, f(x))
	}
}

// Автопоиск интервалов смены знака
func findRootIntervals(a, b float64, m int) [][2]float64 {
	h := (b - a) / float64(m)
	var intervals [][2]float64
	prevX := a
	prevY := f(a)
	for i := 1; i <= m; i++ {
		x := a + float64(i)*h
		y := f(x)
		if prevY*y < 0 {
			intervals = append(intervals, [2]float64{prevX, x})
		}
		prevX = x
		prevY = y
	}
	return intervals
}

// Чтение ввода пользователя
func readFloat(prompt string) float64 {
	reader := bufio.NewReader(os.Stdin)
	for {
		fmt.Print(prompt)
		line, _ := reader.ReadString('\n')
		line = strings.TrimSpace(line)
		v, err := strconv.ParseFloat(line, 64)
		if err == nil {
			return v
		}
		fmt.Println("  Ошибка: введите число.")
	}
}

func main() {
	fmt.Println("=====================================================================")
	fmt.Println("Поиск всех корней функции f(x) = sin^2(x) - x/5 - 1")
	fmt.Println("Метод Ньютона (MN), интервал [-4, 0], точность eps = 1e-4")
	fmt.Println("=====================================================================")

	// 1. Вывод таблицы значений
	printTable(A, B, M)

	// 2. Автопоиск интервалов смены знака (для подсказки)
	intervals := findRootIntervals(A, B, M)
	fmt.Printf("\nНайдено интервалов смены знака: %d\n", len(intervals))
	for i, iv := range intervals {
		fmt.Printf("  %d: [%.4f, %.4f]  f(%.4f)=%.6f, f(%.4f)=%.6f\n",
			i+1, iv[0], iv[1], iv[0], f(iv[0]), iv[1], f(iv[1]))
	}

	// 3. Поиск корней методом Ньютона
	var results []RootResult

	// Если запущено в интерактивном режиме — запрашиваем x0 у пользователя,
	// иначе используем середины интервалов (авто-режим).
	args := os.Args
	interactive := len(args) < 2 || args[1] != "-auto"

	if interactive {
		fmt.Println("\n=== Интерактивный режим ===")
		fmt.Println("Введите начальное приближение x0 для каждого корня (Enter — авто).")
		fmt.Println("Для авто-режима запустите: go run main.go -auto")

		nRoots := 3
		for i := 0; i < nRoots; i++ {
			var x0 float64
			if i < len(intervals) {
				mid := (intervals[i][0] + intervals[i][1]) / 2
				x0 = readFloat(fmt.Sprintf("\nКорень %d, x0 (Enter для авто=%.4f): ", i+1, mid))
				if x0 == 0 && mid != 0 {
					x0 = mid
				}
			} else {
				x0 = readFloat(fmt.Sprintf("\nКорень %d, x0: ", i+1))
			}
			res := newton(x0, Eps, 100)
			results = append(results, res)
			printRootResult(res)
		}
	} else {
		fmt.Println("\n=== Авто-режим ===")
		for i, iv := range intervals {
			x0 := (iv[0] + iv[1]) / 2
			fmt.Printf("\n--- Корень %d, x0 = %.4f ---\n", i+1, x0)
			res := newton(x0, Eps, 100)
			results = append(results, res)
			printRootResult(res)
		}
	}

	// 4. Сохранение данных для графиков
	jsonData, err := json.MarshalIndent(results, "", "  ")
	if err != nil {
		fmt.Println("Ошибка сериализации:", err)
		return
	}
	err = os.WriteFile("roots_data.json", jsonData, 0644)
	if err != nil {
		fmt.Println("Ошибка записи:", err)
	} else {
		fmt.Println("\nДанные сохранены в roots_data.json")
	}
}

func printRootResult(res RootResult) {
	fmt.Printf("\n%-6s | %-14s | %-14s | %-14s\n", "k", "x_k", "f(x_k)", "d_k = x_k - x_{k-1}")
	fmt.Println("------------------------------------------------------------------")
	for _, it := range res.History {
		fmt.Printf("%-6d | %-14.8f | %-14.8e | %-14.8e\n", it.K, it.Xk, it.Fxk, it.Dx)
	}
	fmt.Printf("\n>> Корень: x* = %.8f, f(x*) = %.3e, итераций: %d\n",
		res.Root, res.Froot, res.Iterations)
}
