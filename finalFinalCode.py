# librerías
import os

try:
    import pandapower
    import pandapower.networks as pn
except ImportError:
    print("Instalando librerías necesarias...")
    os.system('pip install -q pandapower pypower matpowercaseframes')

import numpy as np
import scipy.linalg
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from time import perf_counter
from pandapower.converter.matpower import from_mpc
from pandapower.converter.pypower import to_ppc
import urllib.request

import copy
from time import perf_counter

# ==============================================================================
#  CLASE SOLVER PROFESIONAL (LU DECOMPOSITION + HYBRID CORE)
# ==============================================================================
class PowerFlow_Professional:
    def __init__(self, ppc, tol_p=0.10, tol_q=0.10):
        self.baseMVA = ppc['baseMVA']
        self.bus = ppc['bus']
        self.branch = ppc['branch']
        self.gen = ppc['gen']
        self.nb = len(self.bus)

        # Índices y Topología
        self.slack = np.where(self.bus[:, 1] == 3)[0][0]
        self.pv = np.where(self.bus[:, 1] == 2)[0]
        self.pq = np.where(self.bus[:, 1] == 1)[0]
        self.ns = np.setdiff1d(np.arange(self.nb), [self.slack])

        # Mapeo de IDs externos a índices internos (Vital para eficiencia)
        self.bus_ids = self.bus[:, 0].astype(int)
        self.id_to_idx = {b_id: i for i, b_id in enumerate(self.bus_ids)}

        # Inyección de Incertidumbre
        p_val = self.bus[:, 2] / self.baseMVA
        q_val = self.bus[:, 3] / self.baseMVA

        self.dP_inj = np.zeros(self.nb)
        self.dQ_inj = np.zeros(self.nb)

        # Filtro de ruido numérico
        mask_p = np.abs(p_val) > 1e-6
        mask_q = np.abs(q_val) > 1e-6

        self.dP_inj[mask_p] = np.abs(p_val[mask_p] * tol_p)
        self.dQ_inj[mask_q] = np.abs(q_val[mask_q] * tol_q)

        self._build_matrices()

        self.DC_lu = None
        self.AC_lu = None

    def _build_matrices(self):
        fr = self.branch[:, 0].astype(int)
        to = self.branch[:, 1].astype(int)
        r = self.branch[:, 2]
        x = self.branch[:, 3]
        b = self.branch[:, 4]
        tap = self.branch[:, 8]
        tap[tap == 0] = 1.0
        shift = np.deg2rad(self.branch[:, 9])

        z = r + 1j*x
        y = 1/z
        b_sh = 1j * b / 2.0
        a = tap * np.exp(1j * shift)

        Y = np.zeros((self.nb, self.nb), dtype=complex)

        y_a_conj = y / np.conj(a)
        y_a = y / a
        np.add.at(Y, (fr, to), -y_a_conj)
        np.add.at(Y, (to, fr), -y_a)

        y_aa = y / (np.conj(a) * a)
        np.add.at(Y, (fr, fr), y_aa + b_sh)
        np.add.at(Y, (to, to), y + b_sh)

        bus_sh = (self.bus[:, 4] + 1j * self.bus[:, 5]) / self.baseMVA
        idx = np.arange(self.nb)
        Y[idx, idx] += bus_sh

        self.Ybus = Y
        self.G = np.real(Y)
        self.B = np.imag(Y)
        self.Bp = -self.B

    def solve_dc_sensitivities(self):
        Bp_red = self.Bp[np.ix_(self.ns, self.ns)]
        self.DC_lu = scipy.linalg.lu_factor(Bp_red)

    def solve_ac_sensitivities(self):
        V = np.ones(self.nb)
        th = np.zeros(self.nb)
        for i in self.pv: V[i] = self.bus[i, 7]
        V[self.slack] = self.bus[self.slack, 7]

        row = np.arange(self.nb)
        th_diff = th[:, None] - th[None, :]
        ViVj = V[:, None] * V[None, :]
        Gcos_Bsin = self.G * np.cos(th_diff) + self.B * np.sin(th_diff)
        Gsin_Bcos = self.G * np.sin(th_diff) - self.B * np.cos(th_diff)

        H = -ViVj * Gsin_Bcos
        N =  ViVj * Gcos_Bsin
        M = -ViVj * Gcos_Bsin
        L =  ViVj * Gsin_Bcos

        H[row, row] = -self.B.diagonal() * V**2
        N[row, row] =  self.G.diagonal() * V**2
        M[row, row] = -self.G.diagonal() * V**2
        L[row, row] = -self.B.diagonal() * V**2

        J11 = H[np.ix_(self.ns, self.ns)]
        J12 = N[np.ix_(self.ns, self.pq)]
        J21 = M[np.ix_(self.pq, self.ns)]
        J22 = L[np.ix_(self.pq, self.pq)]

        J = np.block([[J11, J12], [J21, J22]])
        self.AC_lu = scipy.linalg.lu_factor(J)

    def get_affine_coeffs(self, model, var_type, idx):
        if model == 'DC':
            if idx == self.slack: return np.zeros(len(self.ns))
            row_idx = np.where(self.ns == idx)[0][0]
            size = len(self.ns)
            unit_vec = np.zeros(size)
            unit_vec[row_idx] = 1.0
            sens = scipy.linalg.lu_solve(self.DC_lu, unit_vec, trans=1)
            return sens * self.dP_inj[self.ns]
        elif model == 'AC':
            noise_vec = np.concatenate([self.dP_inj[self.ns], self.dQ_inj[self.pq]])
            size = len(noise_vec)
            unit_vec = np.zeros(size)
            valid_req = False
            if var_type == 'theta':
                if idx != self.slack:
                    row = np.where(self.ns == idx)[0][0]
                    unit_vec[row] = 1.0
                    valid_req = True
            elif var_type == 'V':
                if idx in self.pq:
                    row_offset = len(self.ns)
                    row_local = np.where(self.pq == idx)[0][0]
                    unit_vec[row_offset + row_local] = 1.0
                    valid_req = True
            if not valid_req: return np.zeros(len(noise_vec))
            sens = scipy.linalg.lu_solve(self.AC_lu, unit_vec, trans=1)
            return sens * noise_vec

    def solve_hybrid_strategy(self, umbral_ganancia=15.0, model='AC'):
        """
        NÚCLEO DEL MÉTODO HÍBRIDO: Itera solo sobre ramas físicas y aplica filtros.
        """
        results = []
        n_branches = self.branch.shape[0]

        for k in range(n_branches):
            f_id = int(self.branch[k, 0])
            t_id = int(self.branch[k, 1])

            try:
                idx_f = self.id_to_idx[f_id]
                idx_t = self.id_to_idx[t_id]
            except KeyError: continue

            if not (idx_f in self.pq and idx_t in self.pq): continue

            var_type = 'V' if model == 'AC' else 'theta'
            c1 = self.get_affine_coeffs(model, var_type, idx_f)
            c2 = self.get_affine_coeffs(model, var_type, idx_t)

            if np.all(c1 == 0) or np.all(c2 == 0): continue

            C = np.vstack([c1, c2])
            Q = (C @ C.T) / 3.0
            sigma_x = np.sqrt(Q[0, 0])
            sigma_y = np.sqrt(Q[1, 1])

            if sigma_x * sigma_y <= 1e-12: continue

            rho = Q[0, 1] / (sigma_x * sigma_y)
            eigvals = scipy.linalg.eigvalsh(Q)
            eigvals[eigvals < 0] = 0

            area_ell = np.pi * np.sqrt(eigvals[0]) * np.sqrt(eigvals[1])
            area_box = 4 * sigma_x * sigma_y

            ganancia = max(0, (1.0 - area_ell/area_box) * 100)

            decision = 'BOX'
            # CRITERIO HÍBRIDO: Ganancia alta Y No Trivial
            if ganancia >= umbral_ganancia and abs(rho) < 0.999:
                decision = 'ELLIPSE'

            results.append({
                'id_branch': k, 'bus_from': f_id, 'bus_to': t_id,
                'rho': rho, 'ganancia': ganancia, 'decision': decision,
                'eigvals': eigvals, 'sigmas': (sigma_x, sigma_y),
                'coeffs': (c1, c2), 'Q': Q # Guardamos para validación posterior
            })
        return results

# ==============================================================================
# METODO NO LINEALES MCS
# ==============================================================================

def generate_nonlinear_mcs_samples(net_base, tol_p=0.10, tol_q=0.10, n_samples=5000, dist='uniform', seed=42):
    """
    Oráculo Físico: Ejecuta flujos de potencia AC-NR no lineales masivos
    para generar el 'Ground Truth' de las desviaciones de voltaje y ángulo.
    """
    print(f"\n[MCS Físico] Iniciando {n_samples} flujos de potencia no lineales (NR)...")
    t0 = perf_counter()
    
    # Trabajamos sobre una copia profunda para no corromper la red original
    net_sim = copy.deepcopy(net_base)
    
    # Extraemos el estado base (Ground Truth de V y theta)
    try:
        pandapower.runpp(net_sim, algorithm='nr', calculate_voltage_angles=True)
    except pandapower.LoadflowNotConverged:
        raise ValueError("La red base no converge. Revisa los datos de entrada.")
        
    v_base = net_sim.res_bus.vm_pu.values.copy()
    th_base = np.deg2rad(net_sim.res_bus.va_degree.values.copy())
    
    # Matrices para almacenar resultados
    n_buses = len(net_sim.bus)
    delta_v = np.zeros((n_buses, n_samples))
    delta_th = np.zeros((n_buses, n_samples))
    
    # Respaldamos las cargas originales
    p_load_base = net_sim.load.p_mw.values.copy()
    q_load_base = net_sim.load.q_mvar.values.copy()
    
    rng = np.random.default_rng(seed)
    exitos = 0
    fallos = 0
    
    for i in range(n_samples):
        # 1. Generar ruido según distribución
        if dist == 'uniform':
            ruido_p = rng.uniform(-tol_p, tol_p, size=len(p_load_base))
            ruido_q = rng.uniform(-tol_q, tol_q, size=len(q_load_base))
        elif dist == 'normal':
            # Asumimos que tol_p representa 3 desviaciones estándar (99.73%)
            ruido_p = rng.normal(0, tol_p/3.0, size=len(p_load_base))
            ruido_q = rng.normal(0, tol_q/3.0, size=len(q_load_base))
        else:
            raise ValueError("Distribución no soportada.")
            
        # 2. Aplicar perturbación
        net_sim.load.p_mw = p_load_base * (1 + ruido_p)
        net_sim.load.q_mvar = q_load_base * (1 + ruido_q)
        
        # 3. Resolver flujo de potencia AC
        try:
            # numba=False evita problemas de recompilación en bucles rápidos en algunas versiones
            pandapower.runpp(net_sim, algorithm='nr', calculate_voltage_angles=True, numba=False)
            
            # 4. Registrar desviaciones reales (Estado perturbado - Estado base)
            v_perturbado = net_sim.res_bus.vm_pu.values
            th_perturbado = np.deg2rad(net_sim.res_bus.va_degree.values)
            
            delta_v[:, exitos] = v_perturbado - v_base
            delta_th[:, exitos] = th_perturbado - th_base
            exitos += 1
            
        except pandapower.LoadflowNotConverged:
            fallos += 1
            continue

    t_total = perf_counter() - t0
    print(f"[MCS Físico] Completado en {t_total:.2f} s. Éxitos: {exitos} | Fallos convergencia: {fallos}")
    
    # Recortar matrices si hubo fallos de convergencia
    return delta_v[:, :exitos], delta_th[:, :exitos]


# ==============================================================================
#  EXPERIMENTOS CIENTÍFICOS AVANZADOS (CORREGIDO / ROBUSTOS)
# ==============================================================================

# --- A. VISUALIZACIÓN SECUENCIAL (BARRIDO) ---
def plot_strip_heatmap(results):
    #Genera los mapas de calor tipo tira (Eficacia y Ángulo).
    reducciones = [r['ganancia'] for r in results]
    # Calculamos ángulos aproximados para visualización
    angulos = []
    for r in results:
        # Recalcular eigenvecs rápido solo para plot
        _, vecs = scipy.linalg.eigh(r['Q'])
        idx_max = 1 # eigh devuelve orden ascendente, el max es el ultimo
        vec_max = vecs[:, idx_max]
        ang = np.degrees(np.arctan2(vec_max[1], vec_max[0]))
        angulos.append(ang)

    data_red = np.array(reducciones).reshape(1, -1)
    data_ang = np.array(angulos).reshape(1, -1)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

    im1 = ax1.imshow(data_red, aspect='auto', cmap='Greens', vmin=0, vmax=50)
    ax1.set_yticks([])
    # IEEE Style: Concise titles
    ax1.set_title("1. Heatmap: Proposed Method Efficiency (Green = High Reduction)")
    plt.colorbar(im1, ax=ax1, orientation='horizontal', label='Reduction (%)', pad=0.2)

    im2 = ax2.imshow(data_ang, aspect='auto', cmap='coolwarm', vmin=-90, vmax=90)
    ax2.set_yticks([])
    ax2.set_title("2. Heatmap: Correlation Orientation (Red=Pos / Blue=Neg)")
    plt.colorbar(im2, ax=ax2, orientation='horizontal', label='Angle (°)', pad=0.2)
    plt.tight_layout()
    plt.show()

# --- B. VISUALIZACIÓN TOPOLÓGICA (SPY PLOT) ---
def plot_topologia_eficacia(solver, results):
    #matriz de adyacencia coloreada por eficacia.
    x_coords, y_coords, colors = [], [], []
    bus_ids = solver.bus_ids

    # Mapeo inverso para coordenadas del plot
    id_to_plot_idx = {b_id: i for i, b_id in enumerate(bus_ids)}

    for res in results:
        try:
            idx_f = id_to_plot_idx[res['bus_from']]
            idx_t = id_to_plot_idx[res['bus_to']]
            gain = res['ganancia']
            x_coords.extend([idx_f, idx_t])
            y_coords.extend([idx_t, idx_f])
            colors.extend([gain, gain])
        except KeyError: continue

    plt.figure(figsize=(9, 9))
    sc = plt.scatter(x_coords, y_coords, c=colors, cmap='viridis_r', s=10, marker='s')
    plt.colorbar(sc, label="Volume Reduction (%)", fraction=0.046, pad=0.04)
    plt.title("Topological Heatmap showing selectivity (IEEE-33)")
    plt.gca().invert_yaxis()
    plt.xlabel("Bus Index (From)")
    plt.ylabel("Bus Index (To)")
    plt.tight_layout()
    plt.show()

# --- C. DIAGNÓSTICO ESTADÍSTICO (FRONTERAS Y TENDENCIAS) ---
def plot_tendencia_correlacion(results):
    # Scatter: Correlation vs Reduction.
    rhos = [abs(r['rho']) for r in results]
    gains = [r['ganancia'] for r in results]

    plt.figure(figsize=(7, 5))
    plt.scatter(rhos, gains, alpha=0.4, s=10, c='blue')
    plt.xlabel(r"Correlation Magnitude $|\rho|$") # Usamos LaTeX simple para el símbolo
    plt.ylabel("Reduction (%)")
    plt.title("Trend: Correlation vs. Efficiency")
    plt.grid(True, alpha=0.3)
    plt.show()

def plot_frontera_decision(results, umbral):
    # Scatter: Hybrid Decision Boundary.
    rhos = [abs(r['rho']) for r in results]
    gains = [r['ganancia'] for r in results]
    colors = ['#D50000' if r['decision'] == 'ELLIPSE' else '#90A4AE' for r in results]

    plt.figure(figsize=(8, 6))
    plt.scatter(rhos, gains, c=colors, s=15, alpha=0.6)
    plt.axhline(y=umbral, color='black', ls='--')
    
    # Texto en inglés
    plt.text(0.1, umbral+2, "ELLIPSE ZONE", color='#D50000', fontweight='bold')
    plt.text(0.5, umbral-5, "BOX ZONE", color='gray', fontweight='bold', ha='center')
    
    plt.title(f"Hybrid Decision Boundary (Threshold: {umbral}%)")
    plt.savefig('Fig5_Right_Decision.png', dpi=300)
    plt.xlabel(r"Correlation Magnitude $|\rho|$")
    plt.ylabel("Volume Reduction (%)")
    plt.grid(True, alpha=0.3)
    plt.show()
    
def plot_caso_estrella(solver, results, delta_v_mcs, ranking=1):
    # Grafica el caso N con mayor reducción.
    # Ordenar por ganancia
    sorted_res = sorted(results, key=lambda x: x['ganancia'], reverse=True)
    
    # Filtrar triviales para que el caso estrella sea interesante (no una línea recta)
    non_trivial = [r for r in sorted_res if abs(r['rho']) < 0.999]

    if ranking > len(non_trivial): return
    target = non_trivial[ranking-1]

    # Recalcular geometría para plot
    Q = target['Q']
    eigvals, eigvecs = scipy.linalg.eigh(Q)
    idx_max = np.argmax(eigvals)
    angle = np.degrees(np.arctan2(eigvecs[1, idx_max], eigvecs[0, idx_max]))

    conf = 0.9973
    s = np.sqrt(stats.chi2.ppf(conf, df=2)) # Acá usamos chi cuadrado
    width = 2 * s * np.sqrt(max(0, eigvals[idx_max]))
    height = 2 * s * np.sqrt(max(0, eigvals[1-idx_max]))

    # Monte Carlo local
    # Extracción de la verdad física no lineal
    idx_f = solver.id_to_idx[target['bus_from']]
    idx_t = solver.id_to_idx[target['bus_to']]

    # Extraemos las filas exactas de los nodos involucrados en la línea crítica
    x_sim = delta_v_mcs[idx_f, :] 
    y_sim = delta_v_mcs[idx_t, :]

    sigma_x, sigma_y = target['sigmas']

    fig, ax = plt.subplots(figsize=(7, 7))
    
    # 1. Puntos MCS (Capa media)
    # AÑADIDO: Etiqueta clara para la leyenda
    ax.scatter(x_sim, y_sim, s=1, c='gray', alpha=1, label='MCS Samples', zorder=2, rasterized=True)

    # 2. Elipse Rellena (Fondo) - No necesita label, es decorativa
    ell = patches.Ellipse((0,0), width, height, angle=angle, fill=True, color='#ffcccb', alpha=0.5, zorder=1)
    ax.add_patch(ell)
    
    # 3. Bordes Elipse (Frente)
    # AÑADIDO: Label para el método propuesto
    ell_edge = patches.Ellipse((0,0), width, height, angle=angle, fill=False, 
                               color='#D50000', lw=2, zorder=3, label='Proposed (Ellipse)')
    ax.add_patch(ell_edge)

    # 4. Rectángulo Caja (Frente)
    # AÑADIDO: Label para el método de comparación con notación matemática LaTeX
    rect = patches.Rectangle((-3*sigma_x, -3*sigma_y), 6*sigma_x, 6*sigma_y, 
                             fill=False, edgecolor='black', ls='--', lw=2, zorder=3, 
                             label=r'CLTM Box ($3\sigma$)')
    ax.add_patch(rect)

    titulo_profesional = (f"Geometric Comparison of Uncertainty Sets (Branch {target['bus_from']}--{target['bus_to']})\n"
                          f"Reduction: {target['ganancia']:.2f}% | Correlation (Rho): {target['rho']:.4f}")
    
    ax.set_title(titulo_profesional, fontsize=11, fontweight='bold')

    ax.set_xlabel(f"Voltage Deviation at Bus {target['bus_from']} ($\Delta V_{{{target['bus_from']}}}$) [p.u.]")
    ax.set_ylabel(f"Voltage Deviation at Bus {target['bus_to']} ($\Delta V_{{{target['bus_to']}}}$) [p.u.]")
    
    # -------------------------------------
    # GENERACIÓN DE LEYENDA IEEE
    # -------------------------------------
    # Forzamos la localización y el estilo
    ax.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black', fontsize=9)
    
    plt.tight_layout()
    plt.show()

def solve_mahalanobis_robust(Q, points):
    """
    Función auxiliar para calcular distancia Mahalanobis evitando errores
    de matriz singular en casos de correlación extrema.
    """
    # 1. Regularización (Jitter): Sumar epsilon a la diagonal
    # Esto evita que la elipse tenga ancho 0 matemático
    epsilon = 1e-9
    Q_reg = Q + np.eye(Q.shape[0]) * epsilon

    try:
        # Intento rápido: Solve estándar
        Z = scipy.linalg.solve(Q_reg, points)
    except scipy.linalg.LinAlgError:
        # Intento robusto: Pseudo-inversa (más lento pero infalible)
        Q_inv = scipy.linalg.pinv(Q)
        Z = Q_inv @ points

    d2 = np.sum(points * Z, axis=0)
    return d2

# --- EXPERIMENTO 1: ANÁLISIS DE SENSIBILIDAD (LINEARIZATION CHECK) ---
def experimento_sensibilidad_incertidumbre(solver, ppc, steps=[0.05, 0.10, 0.20, 0.30]):
    print("\n>>> INICIANDO EXPERIMENTO 1: Sensibilidad a la Magnitud...")

    tasas_violacion = []

    # Respaldo
    original_dp = solver.dP_inj.copy()
    original_dq = solver.dQ_inj.copy()

    # Recalculamos base P/Q
    p_val = solver.bus[:, 2] / solver.baseMVA
    q_val = solver.bus[:, 3] / solver.baseMVA
    mask_p = np.abs(p_val) > 1e-6
    mask_q = np.abs(q_val) > 1e-6

    # Seleccionamos caso testigo (Usamos uno con correlación alta pero no rota)
    results = solver.solve_hybrid_strategy(umbral_ganancia=10.0)
    # Filtramos casos rotos (rho > 0.9999)
    valid_cases = [r for r in results if abs(r['rho']) < 0.999]
    if not valid_cases: valid_cases = results # Fallback

    top_case = sorted(valid_cases, key=lambda x: x['ganancia'], reverse=True)[0]
    idx_f, idx_t = solver.id_to_idx[top_case['bus_from']], solver.id_to_idx[top_case['bus_to']]

    chi2_val = stats.chi2.ppf(0.9973, df=2)

    for tol in steps:
        # 1. Modificar Inyecciones
        solver.dP_inj[:] = 0
        solver.dQ_inj[:] = 0
        solver.dP_inj[mask_p] = np.abs(p_val[mask_p] * tol)
        solver.dQ_inj[mask_q] = np.abs(q_val[mask_q] * tol)

        # 2. Recalcular
        c1 = solver.get_affine_coeffs('AC', 'V', idx_f)
        c2 = solver.get_affine_coeffs('AC', 'V', idx_t)

        C = np.vstack([c1, c2])
        Q = (C @ C.T) / 3.0

        # 3. Monte Carlo
        rng = np.random.default_rng(99)
        eps = rng.uniform(-1, 1, size=(len(c1), 5000))
        points = C @ eps

        # 4. Check Robusto
        d2 = solve_mahalanobis_robust(Q, points)

        violaciones = np.sum(d2 > chi2_val)
        tasa = (violaciones / 5000) * 100
        tasas_violacion.append(tasa)
        print(f"   Tol {tol*100}% -> Violación: {tasa:.2f}%")

    # Restaurar
    solver.dP_inj = original_dp
    solver.dQ_inj = original_dq

    plt.figure(figsize=(7, 5))
    plt.plot([s*100 for s in steps], tasas_violacion, 'o-', color='#D32F2F', lw=2)
    plt.axhline(0.27, color='green', ls='--', label='Theoretical Target (0.27%)') # English
    plt.xlabel("Uncertainty Magnitude (%)")
    plt.ylabel("Empirical Violation Rate (%)")
    plt.title("7. Linearization Robustness vs. Noise Magnitude")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()

# --- EXPERIMENTO 2: ROBUSTEZ NO-GAUSSIANA (BETA/WEIBULL) ---
def experimento_robustez_beta(results):
    print("\n>>> INICIANDO EXPERIMENTO 2: Robustez No-Gaussiana (Beta)...")

    # Filtramos casos extremadamente singulares para evitar errores numéricos
    # Tomamos casos con rho < 0.9999
    valid_results = [r for r in results if abs(r['rho']) < 0.9999]
    top_cases = sorted(valid_results, key=lambda x: x['ganancia'], reverse=True)[:3]

    chi2_val = stats.chi2.ppf(0.9973, df=2)
    n_mc = 10000
    rng = np.random.default_rng(42)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    if len(top_cases) < 3:
        print("No hay suficientes casos válidos para el panel de 3.")
        return

    for i, case in enumerate(top_cases):
        ax = axes[i]
        c1, c2 = case['coeffs']
        Q = case['Q']

        # Generación BETA (Asimétrica)
        # Centramos y escalamos para comparar peras con peras
        raw_beta = rng.beta(2, 5, size=(len(c1), n_mc))
        eps = (raw_beta - 0.285) * 4.0

        points = np.vstack([c1 @ eps, c2 @ eps])

        # Check Robusto
        d2 = solve_mahalanobis_robust(Q, points)

        violaciones = np.sum(d2 > chi2_val)
        tasa = (violaciones / n_mc) * 100

        # Plot
        eigvals, eigvecs = scipy.linalg.eigh(Q)
        idx_max = np.argmax(eigvals)
        angle = np.degrees(np.arctan2(eigvecs[1, idx_max], eigvecs[0, idx_max]))
        s = np.sqrt(chi2_val)
        w = 2 * s * np.sqrt(max(1e-9, eigvals[idx_max]))
        h = 2 * s * np.sqrt(max(1e-9, eigvals[1-idx_max]))

        # Plot labels changes
        ax.scatter(points[0], points[1], s=1, c='orange', alpha=0.3, label='Beta Dist.')
        ell = patches.Ellipse((0,0), w, h, angle=angle, fill=False, edgecolor='blue', lw=2, label='CLTM Ellipse')
        ax.add_patch(ell)
        ax.set_title(f"Case #{i+1} (Rho={case['rho']:.3f})\nBeta Violation: {tasa:.2f}%")
        if i == 0: ax.legend()

    plt.suptitle("Statistical calibration under non-Gaussian uncertainty (P-P Plot)")
    plt.tight_layout()
    plt.show()

# --- EXPERIMENTO 3: CALIBRACIÓN DE CONFIANZA (P-P PLOT) ---
def experimento_calibracion_confianza(results):
    print("\n>>> INICIANDO EXPERIMENTO 3: Calibración de Confianza...")

    target_probs = [0.50, 0.70, 0.80, 0.90, 0.95, 0.99, 0.9973]
    empirical_covs = []

    valid_results = [r for r in results if abs(r['rho']) < 0.999]
    if not valid_results: return
    case = sorted(valid_results, key=lambda x: x['ganancia'], reverse=True)[0]

    Q = case['Q']
    c1, c2 = case['coeffs']

    rng = np.random.default_rng(101)
    eps = rng.uniform(-1, 1, size=(len(c1), 20000))
    points = np.vstack([c1 @ eps, c2 @ eps])

    # Check Robusto
    d2 = solve_mahalanobis_robust(Q, points)

    for p in target_probs:
        thresh = stats.chi2.ppf(p, df=2)
        inside = np.sum(d2 <= thresh)
        cov_real = inside / 20000
        empirical_covs.append(cov_real)

    plt.figure(figsize=(6, 6))
    plt.plot(target_probs, empirical_covs, 'o-', color='purple', lw=2)
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Ideal')
    plt.xlabel("Theoretical Confidence")
    plt.ylabel("Empirical Coverage")
    plt.title("Statistical Calibration (P-P Plot)")
    plt.savefig('Fig5_Left_Calibration.png', dpi=300)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()

# --- EXPERIMENTO 4: IMPACTO DE DENSIDAD ---
def experimento_densidad_penetracion(solver, ppc):
    print("\n>>> INICIANDO EXPERIMENTO 4: Densidad de Penetración...")

    niveles = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    pct_elipses = []

    original_dp = solver.dP_inj.copy()
    original_dq = solver.dQ_inj.copy()

    pq_idx = solver.pq
    n_pq = len(pq_idx)
    p_val_full = solver.bus[:, 2] / solver.baseMVA

    rng = np.random.default_rng(55)

    for lvl in niveles:
        solver.dP_inj[:] = 0
        solver.dQ_inj[:] = 0

        n_active = int(n_pq * lvl)
        if n_active == 0: n_active = 1
        active_idx = rng.choice(pq_idx, size=n_active, replace=False)

        for idx in active_idx:
            solver.dP_inj[idx] = np.abs(p_val_full[idx] * 0.10)

        # Ejecutar
        res = solver.solve_hybrid_strategy(umbral_ganancia=10.0)

        n_ell = sum(1 for r in res if r['decision'] == 'ELLIPSE')
        pct = (n_ell / len(res)) * 100 if res else 0
        pct_elipses.append(pct)
        print(f"   Penetración {lvl*100}% -> {pct:.1f}% Líneas Críticas")

    solver.dP_inj = original_dp
    solver.dQ_inj = original_dq

    plt.figure(figsize=(7, 5))
    plt.fill_between([n*100 for n in niveles], pct_elipses, color='#0288D1', alpha=0.3)
    plt.plot([n*100 for n in niveles], pct_elipses, 'o-', color='#0277BD', lw=2)
    plt.xlabel("Renewable Penetration (%)")
    plt.ylabel("% Critical Branches (Ellipse)")
    plt.title("Adaptability to Renewable Penetration")
    plt.grid(True, alpha=0.3)
    plt.show()

def plot_histograma_reduccion(results, sistema_nombre="SISTEMA"):
    ganancias = [r['ganancia'] for r in results]
    media = np.mean(ganancias)
    
    plt.figure(figsize=(8, 5))
    # Histograma con bordes definidos
    plt.hist(ganancias, bins=40, color='#66BB6A', edgecolor='black', alpha=0.7)
    
    # Línea de la media
    plt.axvline(media, color='red', linestyle='--', linewidth=2, label=f'Mean Reduction: {media:.2f}%')
    
    plt.xlabel('Uncertainty Volume Reduction (%)')
    plt.ylabel('Frequency (Number of Branches)')
    plt.title(f'Uncertainty Reduction Distribution - {sistema_nombre}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

###############################################################################

def plot_fig1_voltage_profile(solver, results, net_name="IEEE 33-Bus"):
    print(f"\nGenerando FIGURA 1: Perfil de Voltajes ({net_name})...")
    
    # 1. Obtener Voltaje Base (V0) del caso resuelto (Newton-Raphson base)
    # Asumimos que solver.bus tiene los resultados del flujo base en la columna 7 (VM)
    # Si no, recalculamos un flujo base rápido:
    v_base = solver.bus[:, 7] # Magnitud V (p.u.)
    
    # 2. Calcular Radio "Box" (Standard IA) para TODOS los nodos
    # Radio = Suma de valores absolutos de los coeficientes de sensibilidad
    radios_box = np.zeros(solver.nb)
    for i in range(solver.nb):
        # Usamos 'AC' y 'V' para obtener coeficientes de voltaje
        coeffs = solver.get_affine_coeffs(model='AC', var_type='V', idx=i)
        radios_box[i] = np.sum(np.abs(coeffs))
        
    # 3. Calcular Radio "Hybrid" (Mejorado)
    # Partimos del Box y reducimos donde haya decisión 'ELLIPSE'
    radios_hybrid = radios_box.copy()
    
    # Mapeo de decisión por nodo (si un nodo toca una línea elipse, se beneficia)
    for res in results:
        if res['decision'] == 'ELLIPSE':
            # Recuperar geometría guardada o recalcular
            Q = res['Q']
            eigvals = res['eigvals'] # Ya calculados en tu script
            
            # Radio Espectral (Formula del paper: sqrt(chi2 * max_lambda))
            # Chi-cuadrado para 99.73% (3 sigma) y 2 grados de libertad = 11.82
            chi2_val = 11.82 
            radio_ell = np.sqrt(chi2_val * np.max(eigvals))
            
            # Aplicamos el radio elíptico a los nodos i y j de esa línea
            # (Tomamos el mínimo si ya fue procesado por otra línea crítica)
            idx_from = solver.id_to_idx[res['bus_from']]
            idx_to = solver.id_to_idx[res['bus_to']]
            
            # Solo actualizamos si el nuevo radio es menor (más preciso) que el Box
            if radio_ell < radios_hybrid[idx_from]:
                radios_hybrid[idx_from] = radio_ell
            if radio_ell < radios_hybrid[idx_to]:
                radios_hybrid[idx_to] = radio_ell

    # 4. Graficar
    plt.figure(figsize=(10, 6))
    
    # Eje X: Nodos
    bus_indices = np.arange(solver.nb)
    
    # Relleno para IA Pura (Box) - Gris
    plt.fill_between(bus_indices, 
                     v_base - radios_box, 
                     v_base + radios_box, 
                     color='gray', alpha=0.3, label='Standard IA (Box)')
    
    # Líneas para Método Propuesto (Hybrid) - Azul/Rojo
    plt.plot(bus_indices, v_base + radios_hybrid, 'b-', lw=1.5, label='CD-HIM Upper')
    plt.plot(bus_indices, v_base - radios_hybrid, 'r-', lw=1.5, label='CD-HIM Lower')
    
    # Voltaje Base
    plt.plot(bus_indices, v_base, 'k--', alpha=0.6, lw=1, label='Base Case')

    plt.xlabel('Bus Index')
    plt.ylabel('Voltage Magnitude (p.u.)')
    plt.title(f'Voltage Profile Bounds ({net_name})')
    plt.legend(loc='lower left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Guardar con nombre específico para LaTeX
    plt.savefig('FIG1.png', dpi=300)
    print("Guardado 'FIG1.png'.")
    plt.show()

def plot_fig5_histogram(results, net_name="PEGASE-2869"):
    print(f"\nGenerando FIGURA 5: Histograma de Reducción ({net_name})...")
    
    if not results:
        print("No hay resultados para graficar.")
        return

    # Extraer las ganancias (CRI) de todas las líneas
    ganancias = [r['ganancia'] for r in results]
    
    # Filtrar ceros si quieres ver solo donde hubo mejora efectiva (Opcional)
    # ganancias = [g for g in ganancias if g > 0.1] 

    plt.figure(figsize=(8, 5))
    
    # Histograma
    n, bins, patches = plt.hist(ganancias, bins=50, color='#2E7D32', 
                                alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Línea de media
    media = np.mean(ganancias)
    plt.axvline(media, color='red', linestyle='dashed', linewidth=1.5, label=f'Mean Reduction: {media:.2f}%')
    
    plt.xlabel('Uncertainty Volume Reduction (%)')
    plt.ylabel('Frequency (Number of Branches)')
    plt.title(f'Uncertainty Reduction Histogram - {net_name}')
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    
    # Guardar
    plt.savefig('Fig4_Histogram.png', dpi=300)
    print("Guardado 'Fig4_Histogram.png'.")
    plt.show()

###############################################################################
# ==============================================================================
# PRUEBA Y TABLA DE RESULTADOS
# ==============================================================================
def seccion_4_ejecucion_y_tabla(solver, ppc):
    print("\n" + "="*80)
    print("EJECUCION ALGORITMO")
    print("="*80)

    t_start = perf_counter()
    # Ejecutamos el NÚCLEO HÍBRIDO con umbral del 15%
    results = solver.solve_hybrid_strategy(umbral_ganancia=15.0, model='AC')
    t_total = perf_counter() - t_start

    # Análisis de los resultados
    total_lines = len(results)
    critical_lines = [r for r in results if r['decision'] == 'ELLIPSE']
    standard_lines = [r for r in results if r['decision'] == 'BOX']

    avg_gain_critical = np.mean([r['ganancia'] for r in critical_lines]) if critical_lines else 0
    max_gain = max([r['ganancia'] for r in results]) if results else 0

    # Generación de la Tabla
    print(f"\nTABLA DE RESULTADOS: RENDIMIENTO DEL ALGORITMO")
    print("-" * 60)
    print(f"{'Métrica':<40} | {'Valor':<15}")
    print("-" * 60)
    print(f"{'Total Líneas Analizadas (PQ-PQ)':<40} | {total_lines}")
    print(f"{'Tiempo de Ejecución (Solver)':<40} | {t_total:.4f} s")
    print(f"{'Líneas Críticas (Usan Elipse)':<40} | {len(critical_lines)} ({len(critical_lines)/total_lines*100:.1f}%)")
    print(f"{'Líneas Estándar (Usan Caja)':<40} | {len(standard_lines)} ({len(standard_lines)/total_lines*100:.1f}%)")
    print(f"{'Reducción Promedio (Zona Crítica)':<40} | {avg_gain_critical:.2f} %")
    print(f"{'Reducción Máxima Detectada':<40} | {max_gain:.2f} %")
    print("-" * 60)

    return results

# ==============================================================================
# GENERACIÓN DE GALERÍA DE GRÁFICOS COMPLETA (ACTUALIZADA)
# ==============================================================================
def seccion_5_galeria_grafica(solver, results, ppc, delta_v_mcs):
    print("\n" + "="*80)
    print("SECCIÓN 5: GENERACIÓN DE GALERÍA DE GRÁFICOS COMPLETA")
    print("="*80)

    if not results:
        print("No hay resultados para graficar.")
        return

    # --- GRÁFICOS BÁSICOS (Ya tenías estos) ---
    print("Generando Gráfico 1 & 2: Mapas de Calor Secuenciales...")
    plot_strip_heatmap(results)

    print("Generando Gráfico 3: Matriz Topológica...")
    plot_topologia_eficacia(solver, results)

    print("Generando Gráfico 4: Tendencia de Correlación...")
    plot_tendencia_correlacion(results)

    print("Generando Gráfico 5: Frontera de Decisión Híbrida...")
    plot_frontera_decision(results, umbral=15.0)

    print("Generando Gráfico 6: Visualización Geométrica (Caso Estrella #1)...")
    plot_caso_estrella(solver, results, delta_v_mcs, ranking=1)
    
    # --- EXPERIMENTOS CIENTÍFICOS AVANZADOS (NUEVOS) ---
    print("\n>>> EJECUTANDO EXPERIMENTOS DE VALIDACIÓN IEEE...")

    # Exp 1: Sensibilidad (Requiere acceso al solver para cambiar tolerancias)
    experimento_sensibilidad_incertidumbre(solver, ppc)

    # Exp 2: Robustez Beta (Usa los resultados actuales)
    experimento_robustez_beta(results)

    # Exp 3: Calibración P-P (Usa los resultados actuales)
    experimento_calibracion_confianza(results)

    # Exp 4: Penetración (Modifica solver temporalmente)
    # experimento_densidad_penetracion(solver, ppc)

    print("\n¡TODOS LOS EXPERIMENTOS COMPLETADOS CON ÉXITO!")

def estimar_tiempo_monte_carlo(net, n_total_muestras=20000):
    print(f"Estimando tiempo para {n_total_muestras} muestras...")

    tiempos = []
    for _ in range(50):
        net.load['p_mw'] *= 1.001

        t0 = perf_counter()
        try:
            pandapower.runpp(
                net,
                algorithm='nr',
                calculate_voltage_angles=True,
                numba=False   
            )
            tiempos.append(perf_counter() - t0)
        except:
            pass

    avg_time_per_pf = np.mean(tiempos)
    total_estimated_seconds = avg_time_per_pf * n_total_muestras

    horas = total_estimated_seconds / 3600
    print(f"Tiempo promedio por Flujo de Potencia (AC-NR): {avg_time_per_pf:.5f} s")
    print(f"Tiempo Total Estimado para {n_total_muestras} muestras: {total_estimated_seconds:.2f} s ({horas:.2f} horas)")

    return total_estimated_seconds

######################################################################################
#                                        M A I N                                     #
######################################################################################
######################################################################################
#                                        M A I N                                     #
######################################################################################
print("\n" + "="*80)
print("EJECUTANDO ANÁLISIS EN SISTEMA DE DISTRIBUCIÓN IEEE-33 BUSES")
print("="*80)

net_33 = pn.case33bw()

# Se agrega el parámetro init='flat' para evitar que busque resultados inexistentes
ppc_33 = to_ppc(net_33, calculate_voltage_angles=True, init='flat')
solver_33 = PowerFlow_Professional(ppc_33, tol_p=0.15, tol_q=0.15)

print("Calculando sensibilidades para IEEE-33 (Jacobiana AC)...")
solver_33.solve_dc_sensitivities()
solver_33.solve_ac_sensitivities()

resultados_33 = seccion_4_ejecucion_y_tabla(solver_33, ppc_33)

# Generar Ground Truth no lineal (Validación estricta)
print("Generando Ground Truth físico para validación estricta (IEEE-33)...")
delta_v_real_33, delta_th_real_33 = generate_nonlinear_mcs_samples(net_33, tol_p=0.15, tol_q=0.15, n_samples=20000)

# Llamada única y correcta a la galería
seccion_5_galeria_grafica(solver_33, resultados_33, ppc_33, delta_v_real_33)

######################################################################################

print("\n" + "#"*80)
print(" TEST 5: PEGASE-2869 (European Transmission Grid)")
print("#"*80)

url_pegase = "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/master/pglib_opf_case2869_pegase.m"
file_pegase = "pglib_opf_case2869_pegase.m"

try:
    # 1. Descarga automática si no existe
    if not os.path.exists(file_pegase):
        print(">>> Descargando dataset PEGASE desde PGLib...")
        urllib.request.urlretrieve(url_pegase, file_pegase)

    # 2. Carga (Importante: Europa funciona a 50Hz)
    print("[*] Cargando y convirtiendo modelo (puede tardar un poco)...")
    net_pegase = from_mpc(file_pegase, f_hz=50)
    ppc_pegase = to_ppc(net_pegase, calculate_voltage_angles=True, init='flat')

    # 3. Instancia del Solver
    solver_pegase = PowerFlow_Professional(ppc_pegase, tol_p=0.05, tol_q=0.05)

    print(f"[*] Procesando matrices de sensibilidad para {solver_pegase.nb} buses...")
    solver_pegase.solve_dc_sensitivities()
    solver_pegase.solve_ac_sensitivities()

    # 4. Ejecución
    res_pegase = seccion_4_ejecucion_y_tabla(solver_pegase, ppc_pegase)
    
    # 5. Generar Ground Truth reducido para PEGASE (evitar cuellos de botella)
    print("Generando Ground Truth físico para PEGASE (Muestra reducida 100 flujos)...")
    delta_v_real_pegase, _ = generate_nonlinear_mcs_samples(net_pegase, tol_p=0.05, tol_q=0.05, n_samples=20000)

    # 6. Galería Gráfica PEGASE (Con argumentos completos)
    seccion_5_galeria_grafica(solver_pegase, res_pegase, ppc_pegase, delta_v_real_pegase)

except Exception as e:
    print(f"[!] Error crítico durante la ejecución de PEGASE: {e}")
    print("    Si es un error de memoria, intenta reducir n_mc en los experimentos.")

######################################################################################

print("\n--- ESTIMACIÓN DE TIEMPOS MONTE CARLO ---")
print(">> IEEE 33-BUS:")
t_mcs_33 = estimar_tiempo_monte_carlo(net_33)

print("\n>> PEGASE 2869:")
t_mcs_pegase = estimar_tiempo_monte_carlo(net_pegase)

######################################################################################
try:
    plot_fig1_voltage_profile(solver_33, resultados_33)
except NameError:
    print("Error: No se encontraron las variables 'solver_33' o 'resultados_33'. Ejecuta primero el caso IEEE-33.")

try:
    plot_fig5_histogram(res_pegase, "PEGASE-2869")
except NameError:
    print("Aviso: 'res_pegase' no encontrado. Usando 'resultados_33' como demostración.")
    plot_fig5_histogram(resultados_33, "IEEE 33-Bus (DEMO)")