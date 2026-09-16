package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

const (
	A         = 0.0
	B         = 4.0
	NumPoints = 21
)

func f(x float64) float64        { return math.Pow(math.Sin(x), 2) - x/5.0 }
func dfExact(x float64) float64  { return math.Sin(2*x) - 0.2 }
func d2fExact(x float64) float64 { return 2 * math.Cos(2*x) }
func integralExact(x float64) float64 {
	return x/2.0 - math.Sin(2*x)/4.0 - x*x/10.0
}

func DF(y []float64, h float64) []float64 {
	n := len(y)
	res := make([]float64, n)
	res[0] = (-3*y[0] + 4*y[1] - y[2]) / (2 * h)
	for i := 1; i < n-1; i++ {
		res[i] = (y[i+1] - y[i-1]) / (2 * h)
	}
	res[n-1] = (3*y[n-1] - 4*y[n-2] + y[n-3]) / (2 * h)
	return res
}

func D2F(y []float64, h float64) []float64 {
	n := len(y)
	res := make([]float64, n)
	for i := range res {
		res[i] = math.NaN()
	}
	for i := 1; i < n-1; i++ {
		res[i] = (y[i-1] - 2*y[i] + y[i+1]) / (h * h)
	}
	return res
}

func MS(a, b float64, m int, fn func(float64) float64) float64 {
	h := (b - a) / float64(m)
	s := 0.0
	for i := 0; i < m; i++ {
		x1 := a + float64(i)*h
		x2 := a + float64(i+1)*h
		s += (fn(x1) + fn(x2)) / 2.0 * h
	}
	return s
}

// --- Структуры данных ---

type DerivData struct {
	H        float64    `json:"h"`
	X        []float64  `json:"x"`
	FExact   []float64  `json:"f_exact"`
	D1Exact  []float64  `json:"d1_exact"`
	D1Approx []float64  `json:"d1_approx"`
	D1Error  []float64  `json:"d1_error"`
	D2Exact  []float64  `json:"d2_exact"`
	D2Approx []*float64 `json:"d2_approx"` // может быть null
	D2Error  []*float64 `json:"d2_error"`  // может быть null
}

type IntegralData struct {
	M      int     `json:"m"`
	H      float64 `json:"h"`
	Approx float64 `json:"approx"`
	Exact  float64 `json:"exact"`
	Error  float64 `json:"error"`
}

type OutputData struct {
	Derivatives []DerivData    `json:"derivatives"`
	Integrals   []IntegralData `json:"integrals"`
}

// Хелпер: NaN/Inf -> nil, иначе указатель
func ptrOrNil(v float64) *float64 {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return nil
	}
	return &v
}

func main() {
	fmt.Println("=====================================================================")
	fmt.Println("Численное дифференцирование и интегрирование")
	fmt.Println("f(x) = sin^2(x) - x/5, [a,b] = [0, 4]")
	fmt.Println("=====================================================================")

	exactIntegral := integralExact(B) - integralExact(A)
	fmt.Printf("\nТочное значение интеграла: %.8f\n", exactIntegral)

	outData := OutputData{}

	hValues := []float64{0.2, 0.1, 0.05}
	nSteps := []int{20, 40, 80}

	for idx, h := range hValues {
		m := nSteps[idx]
		fmt.Printf("\n--- Расчёт производных при h = %.2f (m = %d) ---\n", h, m)

		xTab := make([]float64, m+1)
		yTab := make([]float64, m+1)
		for i := 0; i <= m; i++ {
			xTab[i] = A + float64(i)*h
			yTab[i] = f(xTab[i])
		}

		d1Approx := DF(yTab, h)
		d2Approx := D2F(yTab, h)

		dd := DerivData{H: h}
		step := (B - A) / 20.0
		fmt.Printf("%-4s | %-10s | %-12s | %-12s | %-12s | %-12s\n",
			"j", "x", "f'(x)_точн", "f'(x)_числ", "f''(x)_точн", "f''(x)_числ")
		fmt.Println("----------------------------------------------------------------------------------")

		for j := 0; j < NumPoints; j++ {
			xj := A + float64(j)*step
			idxTab := int(math.Round((xj - A) / h))
			if idxTab < 0 {
				idxTab = 0
			}
			if idxTab > m {
				idxTab = m
			}

			d1e := dfExact(xj)
			d2e := d2fExact(xj)

			var d1a, d2a float64
			if idxTab >= 0 && idxTab <= m {
				d1a = d1Approx[idxTab]
				d2a = d2Approx[idxTab]
			}

			dd.X = append(dd.X, xj)
			dd.FExact = append(dd.FExact, f(xj))
			dd.D1Exact = append(dd.D1Exact, d1e)
			dd.D1Approx = append(dd.D1Approx, d1a)
			dd.D1Error = append(dd.D1Error, d1e-d1a)
			dd.D2Exact = append(dd.D2Exact, d2e)
			dd.D2Approx = append(dd.D2Approx, ptrOrNil(d2a))
			dd.D2Error = append(dd.D2Error, ptrOrNil(d2e-d2a))

			// Красивый вывод (показываем NaN как "---")
			var d2aStr string
			if math.IsNaN(d2a) {
				d2aStr = "---"
			} else {
				d2aStr = fmt.Sprintf("%-12.6f", d2a)
			}

			fmt.Printf("%-4d | %-10.4f | %-12.6f | %-12.6f | %-12.6f | %-12s\n",
				j+1, xj, d1e, d1a, d2e, d2aStr)
		}

		outData.Derivatives = append(outData.Derivatives, dd)
	}

	fmt.Println("\n--- Расчёт интеграла методом трапеций ---")
	mValues := []int{10, 20, 40}
	fmt.Printf("%-6s | %-8s | %-14s | %-14s | %-14s\n",
		"m", "h", "Приближ.", "Точное", "Погрешность")
	fmt.Println("------------------------------------------------------------------")

	for _, m := range mValues {
		approx := MS(A, B, m, f)
		h := (B - A) / float64(m)
		err := exactIntegral - approx
		fmt.Printf("%-6d | %-8.4f | %-14.8f | %-14.8f | %-14.8e\n",
			m, h, approx, exactIntegral, err)

		outData.Integrals = append(outData.Integrals, IntegralData{
			M: m, H: h, Approx: approx, Exact: exactIntegral, Error: err,
		})
	}

	jsonData, err := json.MarshalIndent(outData, "", "  ")
	if err != nil {
		fmt.Println("Ошибка сериализации:", err)
		return
	}
	err = os.WriteFile("diff_int_data.json", jsonData, 0644)
	if err != nil {
		fmt.Println("Ошибка записи файла:", err)
	} else {
		fmt.Println("\nДанные сохранены в diff_int_data.json")
	}
}
