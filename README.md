# 🌀 Driven Cavity Flow Solver using Chorin Projection (GPU-Accelerated | GFDM + Meshgrid)

This project solves the classic **2D lid-driven cavity flow** problem using a **structured mesh-based solver**. The velocity field is computed using the **Chorin projection scheme**, with both **GPU acceleration via CuPy** and **CPU implementation via NumPy** Results are compared with expected flow patterns and validated against benchmark characteristics.

---

## 1. Problem Definition

### Geometry & Boundary Conditions
- **Domain**: $\Omega = [0,1] \times [0,1]$
- **Velocity BCs**:
  - **Top wall** ($y=1$): 
    ```math 
    u(x,1) = 16x^2(1-x)^2, \quad v(x,1) = 0
    ```
    _(Parabolic horizontal velocity profile for a smooth inlet-like behavior)_

  - **Other walls**: 
    ```math
    u = v = 0 \quad \text{(No-slip)}
    ```
- **Pressure BC**: $\frac{\partial p}{\partial n} = 0$ on all walls

### Initial Conditions
```math
\begin{cases}
u(x,y,0) = 0 \\
v(x,y,0) = 0 \\
p(x,y,0) = 0 
\end{cases}
```

---

## 🧮 2. Numerical Method (Chorin Projection + GFDM)

### 🔢 Governing Equations

The incompressible Navier–Stokes equations:

$$
\begin{aligned}
&\frac{d\vec{x}}{dt} = \vec{v} &\text{(Particle motion)} \\\\
&\frac{D\vec{v}}{Dt} = -\nabla p + \nu \Delta \vec{v} + \vec{g} &\text{(Momentum)} \\\\
&\nabla \cdot \vec{v} = 0 &\text{(Incompressibility)}
\end{aligned}
$$

---
## 2. 🧮 Core Numerical Steps (Chorin Projection + GFDM)

### 2.1 Advection
$$
\vec{x}^{n+1} = \vec{x}^n + \Delta t \vec{v}^n
$$

### 2.2 Viscous Step
$$
\vec{v}^* = \vec{v}^n + \Delta t \left( \nu \Delta \vec{v}^n + \vec{g}^n \right)
$$

### 2.3 Pressure Correction
$$
\Delta p^{n+1} = \frac{\nabla \cdot \vec{v}^*}{\Delta t}
$$

### 2.4 Velocity Projection
→ **vⁿ⁺¹ = v* − Δt ∇pⁿ⁺¹**
---

## 🔧 3. Simulation Parameters

| Parameter             | Value     | Description                        |
|-----------------------|-----------|------------------------------------|
| Re                    | 1000      | Reynolds number (implied via ν)    |
| L, K                  | 1.0       | Domain height and width            |
| Nx × Ny               | 41 × 41   | Grid resolution (structured grid)  |
| Δx, Δy                | 0.025     | Grid spacing (uniform)             |
| Δt                    | 0.0005    | Time step                          |
| Final Time            | 5.1 s     | Total simulation duration          |
| Stencil Size          | 9         | Number of neighbors for GFDM       |
| Tracer Points         | 41²       | For visualization of advection     |
| ν                     | 0.001     | Kinematic viscosity (Re=1000)      |

---

## 🧠 4. Derivatives with GFDM (Least Squares)

This solver uses **Generalized Finite Difference Method (GFDM)** to approximate derivatives without a mesh:

```math
\left[\begin{matrix}
f \\ \frac{∂f}{∂x} \\ \frac{∂f}{∂y} \\ \frac{∂²f}{∂x²} \\ \frac{∂²f}{∂x∂y} \\ \frac{∂²f}{∂y²}
\end{matrix}\right]
≈ (M^T W M)^{-1} M^T W
\cdot
\left[\begin{matrix} f_1 \\ f_2 \\ \vdots \\ f_m \end{matrix}\right]
```

### Weights (Gaussian Kernel)
w<sub>i</sub> = exp(−6.25 · ‖x<sub>i</sub> − x‖<sup>2</sup> / h<sup>2</sup>)
---
## ⚙️ 4. Solver Algorithm (Chorin Projection)

1. **Grid Generation**
   - Structured 2D mesh (41 × 41)
   - Boundary vs interior nodes identified

2. **Precompute GFDM Weights**
   - `w_dx`, `w_dy`, `w_lap` (derivative operators)

3. **Initialize Fields**
   ```math
   \begin{cases}
   u(x,y) = 0 \\
   v(x,y) = 0 \\
   p(x,y) = 0 \\
   u(x,1) = 16x^2(1-x)^2 \quad \text{(top wall)}
   \end{cases}


4. **Time Integration**
    - Intermediate velocity:
    ```math
        \vec{u}^* = \vec{u}^n + \Delta t \left(-\vec{u}^n \cdot \nabla \vec{u}^n + \nu \nabla^2 \vec{u}^n\right)
    ```
    - Pressure solve:
    ```math
        \nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \vec{u}^*
    ```
    - Velocity correction:
    ``` math
        \vec{u}^{n+1} = \vec{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}
    ```
    - Enforce boundary conditions

5. **Particle Advection**
   - Tracers evolve using interpolated velocity fields
   - Snapshots stored at \( t = 1.0 \) and \( t = 5.0 \)

---


## 📈 5. Visualization

| Output Field     | Description                        |
|------------------|------------------------------------|
| Tracer Positions | Scatter plot at selected time steps |
| Velocity         | Used for interpolating tracers     |
| Pressure         | Implicit via correction            |


## ✅ 6. Validation Metrics

### Flow Characteristics
- **Primary vortex**: Position near (0.5, 0.75)
- **Secondary vortices**: Weak recirculation zones in corners
  - Bottom-left corner vortex
  - Bottom-right corner vortex

### Quantitative Measures
1. **Divergence-free condition**:
   ```math
   \max|\nabla \cdot \vec{v}| < 10^{-4}
   ```
2. **Steady-state convergence**:
    ```math
    \frac{\|\vec{v}^{n+1} - \vec{v}^n\|}{\|\vec{v}^n\|} < 10^{-6}
    ```
3. **Benchmark comparison**:
    - Primary vortex center position vs. Ghia et al. (1982)
    - Velocity profiles along centerlines

📊 **Performance**:
    - GPU acceleration (CuPy) resulted in a 6× speedup vs NumPy baseline


## 🚀 6. Project Setup & Usage

### 🔧 Requirements

-  Python 3.11+
-  CUDA-enabled GPU (with CUDA 12.x drivers)
-  Anaconda (recommended)

---

### 🧱 Step 1: Clone the repository

```bash
git clone https://github.com/nsreelekha/poisson_equation.git
cd poisson_equation
```
### 🧪 Step 2A: Create and activate the Conda environment
```bash
conda env create -f environment.yml
conda activate cupy129env
```
### 🧪 Step 2B: (Alternative) Use pip + virtualenv
```bash
python -m venv venv
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

pip install -r requirements.txt
```

## 📦 Dependencies
All dependencies are managed via Conda or requirements.txt
### Key packages include:
- cupy-cuda12x — GPU array library (NumPy-compatible)
- numpy, scipy, matplotlib
- pycuda (optional)
- Python 3.11

Install automatically via:
```bash
conda env create -f environment.yml
```

## 🛠️ Troubleshooting
- No GPU?
Replace:
```bash
import cupy as cp
```
with:
```bash
import numpy as cp 
```

- Environment errors?
Re-create the environment:
```bash
conda env remove -n cupy129env
conda env create -f environment.yml
```
- CuPy errors or crashes?
  Ensure that:
    - CUDA 12.x is correctly installed
    - Your GPU is compatible
    - cupy-cuda12x matches your CUDA version
- Simulation unstable?
    Try reducing Δt or increasing Nx, Ny

### 📜 License
This project is licensed under the MIT License.


## 👤 Author

**Sreelekha Nampally**  
🎓 B.Tech in Mathematics and Computing, NIT Mizoram  
🛠️ Project developed as part of Summer Internship at IIT Tirupati  
🌐 [LinkedIn](https://www.linkedin.com/in/sreelekha-nampally) | [GitHub](https://github.com/nsreelekha)
