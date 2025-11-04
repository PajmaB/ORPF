# fpor_lpip_builder.py
import math
import numpy as np
import sympy as sp
import mysql.connector

# ---------------------------
# 0) Parâmetros / Conexão DB
# ---------------------------
DB_CONFIG = {
    'host': 'localhost',
    'user': 'root',
    'password': '',
    'database': 'sistema_eletrico'
}

# ---------------------------
# 1) Ler tabelas do banco
# ---------------------------
cnx = mysql.connector.connect(**DB_CONFIG)
cursor = cnx.cursor()

cursor.execute("SELECT Linha, Barra_Origem, Barra_Destino, Gkm, Bkm, Bsh, Tap FROM dadoslinha")
rows_linhas = cursor.fetchall()
colunas_linhas = [d[0] for d in cursor.description]

cursor.execute("SELECT * FROM dadosbarra")
rows_barras = cursor.fetchall()
colunas_barras = [d[0] for d in cursor.description]

cursor.close()
cnx.close()

# Converter para dicionários mais fáceis de usar
linhas = []
for r in rows_linhas:
    linha_dict = dict(zip(colunas_linhas, r))
    # garantir tipos
    linha_dict['Linha'] = int(linha_dict['Linha'])
    linha_dict['Barra_Origem'] = int(linha_dict['Barra_Origem'])
    linha_dict['Barra_Destino'] = int(linha_dict['Barra_Destino'])
    linha_dict['Gkm'] = float(linha_dict['Gkm'])
    linha_dict['Bkm'] = float(linha_dict['Bkm'])
    linha_dict['Bsh'] = float(linha_dict.get('Bsh', 0.0))
    linha_dict['Tap'] = float(linha_dict.get('Tap', 1.0)) if linha_dict.get('Tap', None) is not None else 1.0
    linhas.append(linha_dict)

# descobrir nome da coluna que representa o id da barra (insensível a maiúsculas)
barra_id_col = None
for c in colunas_barras:
    if c.lower() in ('barra', 'bar', 'id', 'bus', 'bus_id'):
        barra_id_col = c
        break
if barra_id_col is None:
    # fallback: primeira coluna
    barra_id_col = colunas_barras[0]

barras = []
barras_data = {}
for r in rows_barras:
    d = dict(zip(colunas_barras, r))
    bid = int(d[barra_id_col])
    barras.append(bid)
    # garantir campos básicos (com nomes comuns)
    # tenta mapear nomes mais usados: Tipo, PG, PC, QG, QC, Bsh, Vmax, Vmin, Qmax, Qmin
    def get_field(dct, keys, default=0.0):
        for k in keys:
            if k in dct and dct[k] is not None:
                return dct[k]
        return default

    barras_data[bid] = {
        'id': bid,
        'tipo': int(get_field(d, ['Tipo', 'tipo', 'TYPE'], 0)),
        'pg': float(get_field(d, ['PG', 'Pg', 'pg', 'P_G'], 0.0)),
        'pc': float(get_field(d, ['PC', 'Pc', 'pc', 'P_C'], 0.0)),
        'qg': float(get_field(d, ['QG', 'Qg', 'qg', 'Q_G'], 0.0)),
        'qc': float(get_field(d, ['QC', 'Qc', 'qc', 'Q_C'], 0.0)),
        'bsh': float(get_field(d, ['Bsh', 'bsh', 'B_SH'], 0.0)),
        'vmax': float(get_field(d, ['Vmax','V_max','Vmax_k','Vmax_k'] , 1.1)),
        'vmin': float(get_field(d, ['Vmin','V_min','Vmin_k'], 0.9)),
        'qmax': float(get_field(d, ['Qmax','Q_max','QMAX'], 9999.0)),
        'qmin': float(get_field(d, ['Qmin','Q_min','QMIN'], -9999.0)),
    }

# ordenar listas de barras e linhas por id para consistência
barras = sorted(barras)
linhas = sorted(linhas, key=lambda x: x['Linha'])

nb = len(barras)
nL = len(linhas)

# ---------------------------
# 2) Montar Ybus (complex) -> Ybus_G, Ybus_B
# ---------------------------
Ybus = np.zeros((nb, nb), dtype=complex)

# Map barra id -> índice 0-based
barra_to_idx = {b: i for i, b in enumerate(barras)}

for linha in linhas:
    i = barra_to_idx[linha['Barra_Origem']]
    j = barra_to_idx[linha['Barra_Destino']]
    Gkm = linha['Gkm']
    Bkm = linha['Bkm']
    Bsh = linha.get('Bsh', 0.0)
    Tap = linha.get('Tap', 1.0)
    if Tap == 0:
        Tap = 1.0

    # admitância série
    y_series = complex(Gkm, Bkm)  # série (não dividimos por Tap aqui no build simplificado)
    # se for modelar transformador ideal com tap no enrolamento, ajuste aqui. Mantemos forma básica:
    # off-diagonals
    Ybus[i, j] -= y_series
    Ybus[j, i] -= y_series
    # diagonais
    Ybus[i, i] += y_series + complex(0.0, Bsh/2.0)
    Ybus[j, j] += y_series + complex(0.0, Bsh/2.0)

Ybus_G = np.real(Ybus)
Ybus_B = np.imag(Ybus)

# ---------------------------
# 3) Criar variáveis simbólicas (V, Teta, Qg, Tap)
# ---------------------------
# índices simbólicos baseados na lista 'barras' e linhas
V_syms = sp.symbols(f'V1:{nb+1}')       # V1 ... Vnb
Teta_syms = sp.symbols(f'Teta1:{nb+1}') # Teta1 ... Teta_nb
Qg_syms = sp.symbols(f'Qg1:{nb+1}')     # Qg1 ... Qg_nb
Tap_syms = sp.symbols(f'Tap1:{nL+1}')   # Tap1 ... Tap_nL

# map id -> symbol (1-based to match barra ids)
V = {barras[i]: V_syms[i] for i in range(nb)}
Teta = {barras[i]: Teta_syms[i] for i in range(nb)}
Qg = {barras[i]: Qg_syms[i] for i in range(nb)}
Tap = {linhas[i]['Linha']: Tap_syms[i] for i in range(nL)}

# ---------------------------
# 4) Função objetivo simbólica (perdas ativas nas linhas)
#    f(x) = sum_{linha km} Gkm * (Vk^2 + Vm^2 - 2 Vk Vm cos(theta_k - theta_m))
# ---------------------------
funcao_objetivo_fx = 0
for linha in linhas:
    k_id = linha['Barra_Origem']
    m_id = linha['Barra_Destino']
    Gkm = linha['Gkm']
    Vk = V[k_id]
    Vm = V[m_id]
    thk = Teta[k_id]
    thm = Teta[m_id]
    funcao_objetivo_fx += Gkm * (Vk**2 + Vm**2 - 2 * Vk * Vm * sp.cos(thk - thm))

# ---------------------------
# 5) Montagem das restrições de igualdade g(x)
#    - Para cada barra k != slack: Pg - Pc - sum(fluxos) = 0  (ou forma usada: Pk_injetada - pg + pc = 0)
#    - Para barras PQ (tipo == 0): Qk_injetada - qg + qc = 0
#    Observação: aqui usamos a forma Pk_injetada = sum_m Vk*Vm*(Gkm cos + Bkm sin)
# ---------------------------
restricoes_igualdade_gx = []

# se quiser indicar barra slack, ajuste aqui (por enquanto None -> não ignora)
# ex: ignora = 1
ignora = None

for k in barras:
    if k == ignora:
        continue
    # Pk_injetada
    Pk_injetada = 0
    for m in barras:
        Gkm = float(Ybus_G[barra_to_idx[k], barra_to_idx[m]])
        Bkm = float(Ybus_B[barra_to_idx[k], barra_to_idx[m]])
        Pk_injetada += V[k] * V[m] * (Gkm * sp.cos(Teta[k] - Teta[m]) + Bkm * sp.sin(Teta[k] - Teta[m]))
    barra = barras_data[k]
    # equacao: Pk_injetada - pg + pc = 0  (mesma forma que você usou)
    equacao_Pk = sp.simplify(Pk_injetada - barra['pg'] + barra['pc'])
    restricoes_igualdade_gx.append(equacao_Pk)

# Q constraints for PQ buses
for k in barras:
    barra = barras_data[k]
    if int(barra['tipo']) == 0:  # PQ
        Qk_injetada = 0
        for m in barras:
            Gkm = float(Ybus_G[barra_to_idx[k], barra_to_idx[m]])
            Bkm = float(Ybus_B[barra_to_idx[k], barra_to_idx[m]])
            Qk_injetada += V[k] * V[m] * (Gkm * sp.sin(Teta[k] - Teta[m]) - Bkm * sp.cos(Teta[k] - Teta[m]))
        equacao_Qk = sp.simplify(Qk_injetada - barra['qg'] + barra['qc'])
        restricoes_igualdade_gx.append(equacao_Qk)

num_g = len(restricoes_igualdade_gx)

# ---------------------------
# 6) Montagem das restrições de desigualdade h(x) <= 0
#    Tipicamente: V - Vmax <= 0 ; Vmin - V <= 0 ; Qg - Qmax <=0 ; Qmin - Qg <=0 ; Tap - tap_max <=0 ; tap_min - Tap <=0
# ---------------------------
restricoes_desigualdade_hx = []

# Tensão
for k in barras:
    vmax = float(barras_data[k]['vmax'])
    vmin = float(barras_data[k]['vmin'])
    restricoes_desigualdade_hx.append(V[k] - vmax)   # <= 0
    restricoes_desigualdade_hx.append(vmin - V[k])  # <= 0

# Qg limites
for k in barras:
    qmax = float(barras_data[k]['qmax'])
    qmin = float(barras_data[k]['qmin'])
    restricoes_desigualdade_hx.append(Qg[k] - qmax)  # <= 0
    restricoes_desigualdade_hx.append(qmin - Qg[k])  # <= 0

# Tap limites (se tiver valores nas linhas, tentaremos usar campos Tap_min / Tap_max se existirem)
for idx_line, linha in enumerate(linhas):
    # tenta extrair Tap_min / Tap_max se existirem nas colunas; senão usa o valor atual como fixo (sem desiguald)
    tap_max = linha.get('Tap_max', None)
    tap_min = linha.get('Tap_min', None)
    tap_sym = Tap[linha['Linha']]
    if tap_max is not None:
        restricoes_desigualdade_hx.append(tap_sym - float(tap_max))  # <=0
    if tap_min is not None:
        restricoes_desigualdade_hx.append(float(tap_min) - tap_sym)  # <=0

num_h = len(restricoes_desigualdade_hx)

# ---------------------------
# 7) Variáveis de folga e multiplicadores simbólicos
# ---------------------------
s_syms = sp.symbols(f's1:{num_h+1}')
r_syms = sp.symbols(f'r1:{num_h+1}')
lamb_syms = sp.symbols(f'lambda1:{num_g+1}')
pi_syms = sp.symbols(f'pi1:{num_h+1}')

rho, mu = sp.symbols('rho mu')

# ---------------------------
# 8) Montagem da LPIP simbólica
# LPIP = rho * f(x)
#        - mu * sum( log(s_i) + log(r_i) )
#        + sum(r_i)
#        + sum( lambda_j * g_j )
#        + sum( pi_i * (h_i + s_i - r_i) )
# ---------------------------
LPIP = rho * funcao_objetivo_fx
if num_h > 0:
    LPIP = LPIP - mu * sum(sp.log(s_syms[i]) + sp.log(r_syms[i]) for i in range(num_h)) + sum(r_syms)

# adição das somas com multiplicadores
if num_g > 0:
    LPIP = LPIP + sum(lamb_syms[j] * restricoes_igualdade_gx[j] for j in range(num_g))
if num_h > 0:
    LPIP = LPIP + sum(pi_syms[i] * (restricoes_desigualdade_hx[i] + s_syms[i] - r_syms[i]) for i in range(num_h))

LPIP = sp.simplify(LPIP)

print("\n[OK] LPIP simbólica montada.")
print(f" - número de variáveis de estado: V({nb}), Teta({nb}), Qg({nb}), Tap({nL})")
print(f" - número de restrições de igualdade (g): {num_g}")
print(f" - número de restrições de desigualdade (h): {num_h}")

# ---------------------------
# 9) Preparar listas de todas as variáveis para lambdify
#    Ordem: V_syms, Teta_syms, Qg_syms, Tap_syms, s_syms, r_syms, lamb_syms, pi_syms
# ---------------------------
all_vars_syms = list(V_syms) + list(Teta_syms) + list(Qg_syms) + list(Tap_syms) \
                + list(s_syms) + list(r_syms) + list(lamb_syms) + list(pi_syms)

# --- Gradiente simbólico (lista de derivadas parciais)
grad_LPIP = [sp.diff(LPIP, var) for var in all_vars_syms]

# ---------------------------
# 10) Converter para funções numéricas (lambdify)
# LPIP_numeric(x_vector, mu_val, rho_val) where x_vector is same ordering as all_vars_syms
# grad_LPIP_numeric same signature -> returns numpy array
# ---------------------------
# Note: Lambdify with a single vector argument: we'll wrap the arguments inside a function that
# maps the incoming vector into separate symbols (SymPy lambdify doesn't accept list-of-symbols directly
# in a portable way), but we can lambdify using the tuple(all_vars_syms) as arguments and then build wrappers.

LPIP_func = sp.lambdify(tuple(all_vars_syms) + (mu, rho,), LPIP, 'numpy')
grad_funcs = [sp.lambdify(tuple(all_vars_syms) + (mu, rho,), g, 'numpy') for g in grad_LPIP]

def LPIP_numeric(x_vector, mu_val, rho_val):
    """
    x_vector must be a 1D numpy array with length = len(all_vars_syms)
    Order: V1..Vnb, Teta1..Teta_nb, Qg1..Qg_nb, Tap1..Tap_nL, s1..s_num_h, r1..r_num_h, lambda1..lambda_num_g, pi1..pi_num_h
    """
    if len(x_vector) != len(all_vars_syms):
        raise ValueError(f"Expected x_vector length {len(all_vars_syms)}, got {len(x_vector)}")
    args = tuple(x_vector.tolist()) + (mu_val, rho_val)
    return float(LPIP_func(*args))

def grad_LPIP_numeric(x_vector, mu_val, rho_val):
    if len(x_vector) != len(all_vars_syms):
        raise ValueError(f"Expected x_vector length {len(all_vars_syms)}, got {len(x_vector)}")
    args = tuple(x_vector.tolist()) + (mu_val, rho_val)
    grad_vals = np.array([g(*args) for g in grad_funcs], dtype=float)
    return grad_vals

print("\n[OK] Funções numéricas LPIP_numeric e grad_LPIP_numeric prontas para uso.")
print(f"Total de variáveis na LPIP (incluindo folgas e duais): {len(all_vars_syms)}")

# ---------------------------
# 11) (Opcional) Salvar/Exibir expressões importantes para conferência
# ---------------------------
# Exibir função objetivo simbólica resumida:
print("\n--- Função objetivo simbólica (resumida) ---")
sp.pprint(sp.simplify(funcao_objetivo_fx))

# Exibir primeiras 3 restrições g(x)
print("\n--- Primeiras 3 restrições de igualdade g(x) (se existirem) ---")
for i, g in enumerate(restricoes_igualdade_gx[:3]):
    print(f"\ng[{i}] =")
    sp.pprint(sp.simplify(g))

# Exibir primeira restrição h(x) se existir
if num_h > 0:
    print("\n--- Primeira restrição de desigualdade h[0] ---")
    sp.pprint(restricoes_desigualdade_hx[0])

# FIM
