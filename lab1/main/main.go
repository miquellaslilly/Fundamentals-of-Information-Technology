package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

// Problem parameters
const (
	Q   = -2.25
	D   = 4.0
	Eps = 1e-4
	N   = 5 // System order
)

// Matrix represents a 2D slice
type Matrix [][]float64

// Vector represents a 1D slice
type Vector []float64

// Result for a single iteration count (for the chart)
type IterationResult struct {
	Omega     float64 `json:"omega"`
	IterCount int     `json:"iterations"`
}

// Matrix-vector multiplication function
func multiplyMatrixVector(m Matrix, v Vector) Vector {
	size := len(v)
	res := make(Vector, size)
	for i := 0; i < size; i++ {
		sum := 0.0
		for j := 0; j < size; j++ {
			sum += m[i][j] * v[j]
		}
		res[i] = sum
	}
	return res
}

// Compute residual ||Ax - b||
func calculateResidual(m Matrix, x, b Vector) float64 {
	ax := multiplyMatrixVector(m, x)
	maxDiff := 0.0
	for i := range b {
		diff := math.Abs(ax[i] - b[i])
		if diff > maxDiff {
			maxDiff = diff
		}
	}
	return maxDiff
}

// Simple iteration method with relaxation parameter omega
// x_{k+1} = x_k + omega * (b - A*x_k)
func solveSimpleIteration(a Matrix, b Vector, omega float64, eps float64) (Vector, int) {
	n := len(b)
	x := make(Vector, n) // Initial guess x0 = 0

	// Iterative process
	iter := 0
	maxIter := 10000 // Protection against infinite loop

	for {
		iter++

		// Compute Ax
		ax := multiplyMatrixVector(a, x)

		// Compute new approximation
		// x_new = x + omega * (b - ax)
		newX := make(Vector, n)
		diffNorm := 0.0

		for i := 0; i < n; i++ {
			// Formula: x_i^(k+1) = x_i^(k) + omega * (b_i - sum(a_ij * x_j^(k)))
			delta := omega * (b[i] - ax[i])
			newX[i] = x[i] + delta

			if math.Abs(delta) > diffNorm {
				diffNorm = math.Abs(delta)
			}
		}

		// Copy new vector into current one
		copy(x, newX)

		// Convergence check by residual or by solution change
		// Use residual: ||Ax - b|| < eps
		residual := calculateResidual(a, x, b)
		if residual < eps {
			break
		}

		if iter >= maxIter {
			fmt.Printf("Maximum iterations reached for omega=%.2f\n", omega)
			break
		}
	}

	return x, iter
}

func main() {
	// 1. Build matrix A
	A := make(Matrix, N)
	for i := range A {
		A[i] = make(Vector, N)
	}

	// Fill the matrix according to the problem statement
	for i := 0; i < N; i++ {
		// Main diagonal
		if i == 0 || i == N-1 {
			A[i][i] = Q
		} else {
			A[i][i] = -2
		}

		// Subdiagonal
		if i > 0 {
			A[i][i-1] = 1
		}

		// Superdiagonal
		if i < N-1 {
			A[i][i+1] = 1
		}
	}

	// 2. Build right-hand side vector b
	b := make(Vector, N)
	b[0] = 0
	b[1] = D
	b[2] = D
	b[3] = D
	b[4] = 0

	fmt.Println("Matrix A:")
	for _, row := range A {
		fmt.Println(row)
	}
	fmt.Println("Vector b:", b)
	fmt.Println("----------------------------------------")

	// 3. Study the dependence of iteration count on omega
	// Vary omega from 0.2 to 2.0 with step 0.2
	var results []IterationResult
	minIters := 9999999.0
	bestOmega := 0.0

	fmt.Printf("%-10s | %-10s\n", "Omega", "Iterations")
	fmt.Println("---------------------")

	// Generate omega values: 0.2, 0.4, ..., 2.0
	for w := 0.2; w <= 2.0001; w += 0.2 {
		// Round to avoid floating point errors
		omega := math.Round(w*10) / 10

		sol, iters := solveSimpleIteration(A, b, omega, Eps)

		results = append(results, IterationResult{Omega: omega, IterCount: iters})
		fmt.Printf("%-10.1f | %-10d\n", omega, iters)

		if float64(iters) < minIters {
			minIters = float64(iters)
			bestOmega = omega
			// Print the solution for the best omega
			fmt.Printf("  -> Solution at omega=%.1f: [", omega)
			for i, v := range sol {
				fmt.Printf("%.4f", v)
				if i < len(sol)-1 {
					fmt.Print(", ")
				}
			}
			fmt.Println("]")
		}
	}

	fmt.Println("----------------------------------------")
	fmt.Printf("Minimum number of iterations: %d at omega = %.1f\n", int(minIters), bestOmega)

	// 4. Save chart data to JSON
	jsonData, _ := json.MarshalIndent(results, "", "  ")
	err := os.WriteFile("data.json", jsonData, 0644)
	if err != nil {
		fmt.Println("Error writing data.json file:", err)
	} else {
		fmt.Println("Chart data saved to data.json")
	}
}
