import pandas as pd
import numpy as np
import sympy as sp
import math
import mysql.connector

# Conexão com banco de dados
conexao = mysql.connector.connect(
    host='localhost',
    user='root',
    password='',
    database='sistema_eletrico'
)

funcao_objetivo_fx = 0
restricoes_igualdade_gx = []
restricoes_desigualdade_hx = []

cursor = conexao.cursor()

# Quantidade de linhas e barras
cursor.execute("SELECT * FROM dadoslinha")
registros = cursor.fetchall()
quantDeRegistros = len(registros)

cursor.execute("SELECT * FROM dadosbarra")
registros2 = cursor.fetchall()
quantDeBarras = len(registros2)


# Definindo variáveis simbólicas fora do loop (mantidas apenas para a FO)
GkmR, TkmR, VkR, VmR, ThetaKR, ThetaMR, BkmR, BshR = sp.symbols("Gkm Tkm Vk Vm ThetaK ThetaM Bkm Bsh")

# Expressão simbólica da função objetivo (Perdas Ativas)
funcaoResol = "Gkm*((1/Tkm**2)*Vk**2 + Vm**2 - 2*(1/Tkm)*Vk*Vm*cos(ThetaK - ThetaM))"
expr = sp.sympify(funcaoResol, locals={"cos": sp.cos})


resultadoF = 0
FuncaoObjetivoEscrita = ""

for x in range(1, quantDeRegistros + 1):
    cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (x,))
    resultado = cursor.fetchone()
    if not resultado:
        continue

    Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = resultado

    cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
    Vk = cursor.fetchone()[0]
    cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
    Vm = cursor.fetchone()[0]
    cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
    ThetaK = cursor.fetchone()[0]
    cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
    ThetaM = cursor.fetchone()[0]

    # se Tap inválido
    if Tkm == 0: Tkm = 1

    # termo conforme referência: g*(Vk**2 + Vm**2 - 2*Vk*Vm*cos(...))
    termo_val = Gkm * (Vk**2 + Vm**2 - 2 * Vk * Vm * math.cos(ThetaK - ThetaM))
    resultadoF += termo_val

    termo_str = f"{Gkm}*({Vk}**2 + {Vm}**2 - 2*{Vk}*{Vm}*cos({ThetaK} - {ThetaM}))"
    if x < quantDeRegistros:
        termo_str += " + "
    FuncaoObjetivoEscrita += termo_str

print(f"\nFunção Objetivo:\n{FuncaoObjetivoEscrita}")
print("Resultado da FO: ", float(resultadoF))

# Ignorar a barra slack
ignora = int(input("\nQual é a barra Slack? \n"))

# Mapeando as linhas por barra origem e destino
linhasPorBarra_Origem = {barra: [] for barra in range(1, quantDeBarras + 1)}
linhasPorBarra_Destino = {barra: [] for barra in range(1, quantDeBarras + 1)}

cursor.execute("SELECT Linha, Barra_Origem, Barra_Destino FROM dadoslinha")
todas_linhas = cursor.fetchall()
for linha, origem, destino in todas_linhas:
    linhasPorBarra_Origem[origem].append(linha)
    linhasPorBarra_Destino[destino].append(linha)

# -------------------------------------------------------------------------
# FUNÇÕES DE FLUXO DE POTÊNCIA ATIVA (Pkm) - CORREÇÃO JÁ APLICADA
# -------------------------------------------------------------------------

def Pkm_Origem(g, b, tap, Vk, Vm, ThetaK, ThetaM):
    # ((tap*Vk)**2)*g - (tap*Vk)*Vm*( g*cos(ThetaK-ThetaM) + b*sin(ThetaK-ThetaM) )
    return ((tap * Vk) ** 2) * g - (tap * Vk) * Vm * (g * math.cos(ThetaK - ThetaM) + b * math.sin(ThetaK - ThetaM))

def Pkm_Destino(g, b, tap, Vk_destino, Vm_origem, ThetaK_destino, ThetaM_origem):
    # ((V_m)**2)*g - (tap*V_m)*V_k*( g*cos(...) + b*sin(...) )
    # aqui Vk_destino equivale a V_m, Vm_origem equivale a V_k
    return (Vk_destino ** 2) * g - (tap * Vk_destino) * Vm_origem * (g * math.cos(ThetaK_destino - ThetaM_origem) + b * math.sin(ThetaK_destino - ThetaM_origem))


# Reconstruir restrição P (Pkm − PGk + PCk = 0 -> no ref: Pg - (Pc) - (termos) = 0)
funcRest1Escrita = ""

for i in range(1, quantDeBarras + 1):
    if i == ignora:
        continue

    cursor.execute("SELECT Tipo, PG, QG, Qmin_G, Qmax_G, PC, QC, Bsh FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Pg, Qg, Qmin_G, Qmax_G, Pc, Qc, Bsh_barra = cursor.fetchone()

    potencia_fluxo_i = 0
    termos = []

    # Origem = i
    for linha in linhasPorBarra_Origem[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = cursor.fetchone()
        if Tkm == 0: Tkm = 1

        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        Vk = cursor.fetchone()[0]
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        Vm = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        ThetaK = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        ThetaM = cursor.fetchone()[0]

        potencia_fluxo_linha = Pkm_Origem(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM)
        potencia_fluxo_i += potencia_fluxo_linha

        termo = f"+((({Tkm}*{Vk})**2)*{Gkm} - ({Tkm}*{Vk})*{Vm}*({Gkm}*cos({ThetaK}-{ThetaM}) + ({Bkm})*sin({ThetaK}-{ThetaM})))\n"
        termos.append(termo)

    # Destino = i
    for linha in linhasPorBarra_Destino[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = cursor.fetchone()
        if Tkm == 0: Tkm = 1

        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        Vm_origem = cursor.fetchone()[0]
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        Vk_destino = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        ThetaM_origem = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        ThetaK_destino = cursor.fetchone()[0]

        potencia_fluxo_linha = Pkm_Destino(Gkm, Bkm, Tkm, Vk_destino, Vm_origem, ThetaK_destino, ThetaM_origem)
        potencia_fluxo_i += potencia_fluxo_linha

        termo = f"+(({Vk_destino}**2)*{Gkm} - ({Tkm}*{Vk_destino})*{Vm_origem}*({Gkm}*cos({ThetaK_destino}-{ThetaM_origem}) + ({Bkm})*sin({ThetaK_destino}-{ThetaM_origem})))\n"
        termos.append(termo)

    # Ordem igual ao referência: Pg - Pc - (termos) = 0
    resultadoRest1F = Pg - Pc - potencia_fluxo_i
    funcRest1 = " ".join(termos)
    funcRest1Escrita += f"\nRestrição da Barra {i}:\n{Pg} - ({Pc}) - ({funcRest1}) = {resultadoRest1F}\n ----------------------- X --------------------- \n"

# ============================
# Funções Q (usar tap*V igual referência)
# ============================

def Qkm_Origem(g, b, tap, Vk, Vm, ThetaK, ThetaM, bsh_linha):
    # -((tap*Vk)**2)*(b + bsh_linha) + (tap*Vk)*Vm*(b*cos - g*sin)
    return -((tap * Vk) ** 2) * (b + bsh_linha) + (tap * Vk) * Vm * (b * math.cos(ThetaK - ThetaM) - g * math.sin(ThetaK - ThetaM))

def Qkm_Destino(g, b, tap, Vk_destino, Vm_origem, ThetaK_destino, ThetaM_origem, bsh_linha):
    # -(Vk_destino**2)*(b + bsh_linha) + (tap*Vk_destino)*Vm_origem*(b*cos - g*sin)
    return -(Vk_destino ** 2) * (b + bsh_linha) + (tap * Vk_destino) * Vm_origem * (b * math.cos(ThetaK_destino - ThetaM_origem) - g * math.sin(ThetaK_destino - ThetaM_origem))


# Reconstruir restrição Q (forma do ref: Qg + Bsh*(V_i**2) - Qc - (termos) = 0)
funcRest2Escrita = ""

for i in range(1, quantDeBarras + 1):
    if i == ignora:
        continue

    cursor.execute("SELECT Tipo, PG, QG, Qmin_G, Qmax_G, PC, QC, Bsh FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Pg, Qg, Qmin_G, Qmax_G, Pc, Qc, Bshk = cursor.fetchone()

    # Só monta restrição Q para barras PQ (Tipo == 0) no estilo do referência
    if Tipo != 0:
        continue

    potencia_refluxo_i = 0
    termos = []

    # Origem = i
    for linha in linhasPorBarra_Origem[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino, Bsh FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino, Bsh_linha = cursor.fetchone()
        if Tkm == 0: Tkm = 1

        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        Vk = cursor.fetchone()[0]
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        Vm = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        ThetaK = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        ThetaM = cursor.fetchone()[0]

        potencia_refluxo_linha = Qkm_Origem(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM, Bsh_linha)
        potencia_refluxo_i += potencia_refluxo_linha

        termo = f"-(({Tkm}*{Vk})**2)*({Bkm}+{Bsh_linha}) + ({Tkm}*{Vk})*{Vm}*({Bkm}*cos({ThetaK}-{ThetaM}) - {Gkm}*sin({ThetaK}-{ThetaM}))\n"
        termos.append(termo)

    # Destino = i
    for linha in linhasPorBarra_Destino[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino, Bsh FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino, Bsh_linha = cursor.fetchone()
        if Tkm == 0: Tkm = 1

        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        Vm_origem = cursor.fetchone()[0]
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        Vk_destino = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        ThetaM_origem = cursor.fetchone()[0]
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        ThetaK_destino = cursor.fetchone()[0]

        potencia_refluxo_linha = Qkm_Destino(Gkm, Bkm, Tkm, Vk_destino, Vm_origem, ThetaK_destino, ThetaM_origem, Bsh_linha)
        potencia_refluxo_i += potencia_refluxo_linha

        termo = f"-({Vk_destino}**2)*({Bkm}+{Bsh_linha}) + ({Tkm}*{Vk_destino})*{Vm_origem}*({Bkm}*cos({ThetaK_destino}-{ThetaM_origem}) - {Gkm}*sin({ThetaK_destino}-{ThetaM_origem}))\n"
        termos.append(termo)

    # Qsh do barramento i
    cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (i,))
    Vk_i = cursor.fetchone()[0]
    QSHk = Bshk * (Vk_i ** 2)

    # Forma do referência: Qg + Bsh*(V_i**2) - Qc - (termos) = 0
    resultadoRest2F = Qg + QSHk - Qc - potencia_refluxo_i

    restricao2 = " + ".join(termos)
    funcRest2Escrita += f"\nRestrição da Barra {i}:\n{Qg} + {QSHk} - {Qc} - ({restricao2}) = {resultadoRest2F}\n ----------------------- X --------------------- \n"

funcRest3Escrita = ""
for i in range(1, quantDeBarras + 1):
    cursor.execute("SELECT Tipo, Qmin_G, QG, Qmax_G FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Qmin_G, Qg, Qmax_G = cursor.fetchone()
    if Tipo == 0:
        continue
    termo3 = f"{Qmin_G} <= QG{i} <= {Qmax_G}"
    funcRest3Escrita += f"\nRestrição de Qg da Barra {i}:\n{termo3}\n-----------------------------\n"

funcRest4Escrita = ""
for i in range(1, quantDeBarras + 1):
    cursor.execute("SELECT V FROM dadosbarra WHERE Barra = %s", (i,))
    result = cursor.fetchone()
    V = result[0]
    termo4 = f"{10.0} <= {V} <= {10.0}"
    funcRest4Escrita += f"\nRestrição de Tensão da Barra {i}:\n{termo4}\n-----------------------------\n"

funcRest5Escrita = ""
cursor.execute("SELECT Linha, Barra_Origem, Barra_Destino FROM dadoslinha")
todas_linhas = cursor.fetchall()

for linha_id, barraOrigem, barraDestino in todas_linhas:
    termo5 = f"10 <= T{barraOrigem}_{barraDestino} <= 10"
    funcRest5Escrita += f"\nRestrição de Tap da Linha {linha_id} (Barra {barraOrigem} - Barra {barraDestino}):\n{termo5}\n-----------------------------\n"

# Fechando a conexão
cursor.close()
conexao.close()
    
print("_________________________RESTRIÇAO 1_________________________")
print("_________________________Pkm − PGk + PCk = 0, ∀k ∈ G' ∪ C_________________________")
print(funcRest1Escrita)
print("_________________________RESTRIÇAO 2_________________________")
print("_________________________Qkm − QGk + QCk − QSHk = 0, ∀k ∈ C_________________________")
print(funcRest2Escrita)
print("_________________________RESTRIÇAO 3_________________________")
print(funcRest3Escrita)
print("_________________________RESTRIÇAO 4_________________________")
print(funcRest4Escrita)
print("_________________________RESTRIÇAO 5_________________________")
print(funcRest5Escrita)

#----------------------------------------X Otimizador X-------------------------------------

# fpor_lpip_numeric.py
import math
import numpy as np
import mysql.connector
from scipy.optimize import minimize

# --------------------------
# Parâmetros (ajuste se quiser)
# --------------------------
DB_CONFIG = {
    'host': 'localhost',
    'user': 'root',
    'password': '',
    'database': 'sistema_eletrico'
}

# Parâmetros LPIP iniciais
mu_init = 1.0        # parâmetro de barreira inicial
rho_init = 1.0       # fator de escalonamento da função objetivo
mu_reduction = 0.1   # fator para reduzir mu por iteração externa
outer_iters = 4      # número de iterações externas (barreira)
tol_inner = 1e-6     # tolerância do otimizador interno
min_slack = 1e-6     # limite inferior para s e r (evitar log(0))

# Vmin/Vmax e Tap como especificado
Vmin_global = 10.0
Vmax_global = 10.0
Tap_min_global = 10.0
Tap_max_global = 10.0

# --------------------------
# 1) Ler dados do banco e carregar em memória
# --------------------------
cnx = mysql.connector.connect(**DB_CONFIG)
cursor = cnx.cursor()

# Ler linhas
cursor.execute("SELECT Linha, Barra_Origem, Barra_Destino, Gkm, Bkm, Bsh, Tap FROM dadoslinha")
linhas_rows = cursor.fetchall()
linhas_cols = [d[0] for d in cursor.description]

# Ler barras
cursor.execute("SELECT * FROM dadosbarra")
barras_rows = cursor.fetchall()
barras_cols = [d[0] for d in cursor.description]

cursor.close()
cnx.close()

# Montar data structures
linhas = []
for r in linhas_rows:
    d = dict(zip(linhas_cols, r))
    linhas.append({
        'Linha': int(d['Linha']),
        'orig': int(d['Barra_Origem']),
        'dest': int(d['Barra_Destino']),
        'Gkm': float(d['Gkm']),
        'Bkm': float(d['Bkm']),
        'Bsh': float(d.get('Bsh', 0.0) if d.get('Bsh', None) is not None else 0.0),
        'Tap': float(d.get('Tap', 1.0) if d.get('Tap', None) is not None else 1.0)
    })

# Detectar coluna de id da barra (flexível)
barra_id_col = None
for c in barras_cols:
    if c.lower() in ('barra', 'bar', 'id', 'bus', 'bus_id'):
        barra_id_col = c
        break
if barra_id_col is None:
    barra_id_col = barras_cols[0]

barras = []
barras_data = {}
for r in barras_rows:
    d = dict(zip(barras_cols, r))
    bid = int(d[barra_id_col])
    barras.append(bid)
    barras_data[bid] = {
        'Tipo': int(d.get('Tipo', 0) if d.get('Tipo', None) is not None else 0),
        'PG': float(d.get('PG', 0.0) if d.get('PG', None) is not None else 0.0),
        'PC': float(d.get('PC', 0.0) if d.get('PC', None) is not None else 0.0),
        'QG': float(d.get('QG', 0.0) if d.get('QG', None) is not None else 0.0),
        'QC': float(d.get('QC', 0.0) if d.get('QC', None) is not None else 0.0),
        'Qmin_G': float(d.get('Qmin_G', -9999.0) if d.get('Qmin_G', None) is not None else -9999.0),
        'Qmax_G': float(d.get('Qmax_G', 9999.0) if d.get('Qmax_G', None) is not None else 9999.0),
        # Tensão V e Teta iniciais (se a tabela tiver colunas V e Teta)
        'V_init': float(d.get('V', 1.0) if d.get('V', None) is not None else 1.0),
        'Teta_init': float(d.get('Teta', 0.0) if d.get('Teta', None) is not None else 0.0)
    }

barras = sorted(barras)
nb = len(barras)
nL = len(linhas)

# Map barra id -> índice 0-based
idx_of_barra = {barras[i]: i for i in range(nb)}

# Build adjacency matrices Ybus_G and Ybus_B from lines (full Ybus)
Ybus = np.zeros((nb, nb), dtype=complex)
for ln in linhas:
    i = idx_of_barra[ln['orig']]
    j = idx_of_barra[ln['dest']]
    Gkm = ln['Gkm']
    Bkm = ln['Bkm']
    Bsh = ln.get('Bsh', 0.0)
    # série
    y = complex(Gkm, Bkm)
    Ybus[i, j] -= y
    Ybus[j, i] -= y
    Ybus[i, i] += y + complex(0.0, Bsh/2.0)
    Ybus[j, j] += y + complex(0.0, Bsh/2.0)

Ybus_G = np.real(Ybus)
Ybus_B = np.imag(Ybus)

# Para fazer cálculo de perdas, mantemos lista de linhas com Gkm e conectividade
# já temos em 'linhas' com keys orig,dest,Gkm

# --------------------------
# Solicitar barra slack ao usuário (como antes)
# --------------------------
ignora = int(input("Qual é a barra Slack? (digite o número da barra) \n"))

# --------------------------
# Construir índices das variáveis no vetor x
# Ordem: V(1..nb), Teta(1..nb), Qg(1..nb), Tap(1..nL), s(1..nh), r(1..nh), lambda(1..ng), pi(1..nh)
# Primeiro precisamos determinar ng (número de equações de igualdade) e nh (número desigualdades)
# --------------------------

# Número de restrições de igualdade:
# - P balance por barra exceto slack -> nb - 1
# - Q balance para barras PQ (Tipo == 0) -> count_pq
count_pq = sum(1 for b in barras if barras_data[b]['Tipo'] == 0)
ng = (nb - 1) + count_pq

# Número de desigualdades:
# Para cada barra: 2 para V (V - Vmax <=0 e Vmin - V <=0) => 2*nb
# Para cada barra: 2 para Qg limits (Qg - Qmax <=0 e Qmin - Qg <=0) => 2*nb
# For taps: for each line maybe 2, but as user mandated taps 10..10 => still include them for completeness
nh = 2 * nb + 2 * nb  # V + Qg
# add tap inequalities (2 per line) if you want:
nh += 2 * nL

# Build index map
base = 0
idx = {}
idx['V'] = (base, base + nb); base += nb
idx['Teta'] = (base, base + nb); base += nb
idx['Qg'] = (base, base + nb); base += nb
idx['Tap'] = (base, base + nL); base += nL
idx['s'] = (base, base + nh); base += nh
idx['r'] = (base, base + nh); base += nh
idx['lamb'] = (base, base + ng); base += ng
idx['pi'] = (base, base + nh); base += nh

n_vars = base

# Helper para extrair sub-vetores
def get_slice(x, key):
    a, b = idx[key]
    return x[a:b]

# --------------------------
# Função objetivo (numérica) - perdas iguais ao seu resultadoF
# f(x) = sum_{linhas} Gkm*(Vk^2 + Vm^2 - 2 Vk Vm cos(theta_k - theta_m))
# --------------------------
def f_loss(x):
    Vx = get_slice(x, 'V')
    Thetax = get_slice(x, 'Teta')
    total = 0.0
    for ln in linhas:
        k = ln['orig']
        m = ln['dest']
        i = idx_of_barra[k]
        j = idx_of_barra[m]
        Gkm = ln['Gkm']
        Vk = Vx[i]
        Vm = Vx[j]
        thk = Thetax[i]
        thm = Thetax[j]
        total += Gkm * (Vk**2 + Vm**2 - 2.0 * Vk * Vm * math.cos(thk - thm))
    return total

# --------------------------
# Restrições de igualdade g(x) -> retorna vetor (length ng)
# Ordem: primeiro P balances (barra 1..nb except slack), depois Q balances for PQ buses in ascending barra id
# Pk_injetada = sum_m Vk Vm (Gkm cos + Bkm sin)
# return g = [g1, g2, ...] where each g = Pk_injetada - PG + PC
# --------------------------
def g_vector(x):
    Vx = get_slice(x, 'V')
    Thetax = get_slice(x, 'Teta')
    Qgx = get_slice(x, 'Qg')
    # Build dictionary of P balances
    g_list = []
    # P balances (skip slack)
    for k in barras:
        if k == ignora:
            continue
        Pk = 0.0
        i = idx_of_barra[k]
        for m in barras:
            j = idx_of_barra[m]
            Gkm = Ybus_G[i, j]
            Bkm = Ybus_B[i, j]
            Pk += Vx[i] * Vx[j] * (Gkm * math.cos(Thetax[i] - Thetax[j]) + Bkm * math.sin(Thetax[i] - Thetax[j]))
        Pg = barras_data[k]['PG']
        Pc = barras_data[k]['PC']
        g_val = Pk - Pg + Pc  # note: matches equacao_Pk = Pk_injetada - pg + pc
        g_list.append(g_val)
    # Q balances for PQ buses
    for k in barras:
        if barras_data[k]['Tipo'] == 0:
            Qk = 0.0
            i = idx_of_barra[k]
            for m in barras:
                j = idx_of_barra[m]
                Gkm = Ybus_G[i, j]
                Bkm = Ybus_B[i, j]
                Qk += Vx[i] * Vx[j] * (Gkm * math.sin(Thetax[i] - Thetax[j]) - Bkm * math.cos(Thetax[i] - Thetax[j]))
            Qg = Qgx[i]  # variável Qg
            Qc = barras_data[k]['QC']
            gQ = Qk - Qg + Qc  # note: equacao_Qk = Qk_injetada - qg + qc
            g_list.append(gQ)
    return np.array(g_list, dtype=float)

# --------------------------
# Restrições de desigualdade h(x) <= 0  -> retorna vetor length nh
# Order consistent with nh definition:
# 1..nb: V - Vmax <= 0
# nb+1..2nb: Vmin - V <= 0
# 2nb+1..3nb: Qg - Qmax <=0
# 3nb+1..4nb: Qmin - Qg <=0
# following 2*nL: Tap - Tap_max <=0 and Tap_min - Tap <=0
# --------------------------
def h_vector(x):
    Vx = get_slice(x, 'V')
    Qgx = get_slice(x, 'Qg')
    Tapx = get_slice(x, 'Tap')
    h = []
    # V upper
    for i, b in enumerate(barras):
        h.append(Vx[i] - Vmax_global)
    # V lower
    for i, b in enumerate(barras):
        h.append(Vmin_global - Vx[i])
    # Qg upper
    for i, b in enumerate(barras):
        qmax = barras_data[b]['Qmax_G']
        h.append(Qgx[i] - qmax)
    # Qg lower
    for i, b in enumerate(barras):
        qmin = barras_data[b]['Qmin_G']
        h.append(qmin - Qgx[i])
    # Tap upper/lower per line (use global bounds given)
    for ln in linhas:
        # Tap - Tap_max <= 0
        h.append(Tapx[ln['Linha'] - 1] - Tap_max_global)
    for ln in linhas:
        # Tap_min - Tap <=0
        h.append(Tap_min_global - Tapx[ln['Linha'] - 1])
    return np.array(h, dtype=float)

# --------------------------
# LPIP function numeric
# LPIP = rho * f(x) - mu * sum(log(s)+log(r)) + sum(r) + lamb^T g + pi^T (h + s - r)
# x includes full vector (primal + slacks + duais)
# --------------------------
def LPIP_value(x, mu_val, rho_val):
    # extract pieces
    s = get_slice(x, 's')
    r = get_slice(x, 'r')
    lamb = get_slice(x, 'lamb')
    pi = get_slice(x, 'pi')
    # primal parts are in x as well
    # compute f(x)
    fx = f_loss(x)
    gx = g_vector(x)
    hx = h_vector(x)
    # barrier: ensure s and r > 0 (they are bounded in optimizer but numeric safety)
    if np.any(s <= 0) or np.any(r <= 0):
        # valor alto penaliza
        return 1e50 + np.sum((s <= 0).astype(float)) * 1e6 + np.sum((r <= 0).astype(float)) * 1e6
    term_barreira = - mu_val * (np.sum(np.log(s)) + np.sum(np.log(r)))
    val = rho_val * fx + term_barreira + np.sum(r) + np.dot(lamb, gx) + np.dot(pi, (hx + s - r))
    return float(val)

# --------------------------
# Função wrapper para o otimzador (retorna scalar)
# --------------------------
def objective_for_optimizer(x, mu_val, rho_val):
    return LPIP_value(x, mu_val, rho_val)

# --------------------------
# Montar x0 (chute inicial)
# --------------------------
def build_initial_x():
    x0 = np.zeros(n_vars)
    # V initial: use V_init from table or 1.0 (but user set Vmin/Vmax 10: set to 10)
    for i, b in enumerate(barras):
        x0[idx['V'][0] + i] = 10.0  # conforme solicitado (10)
    # Theta initial
    for i, b in enumerate(barras):
        x0[idx['Teta'][0] + i] = barras_data[b]['Teta_init'] if 'Teta_init' in barras_data[b] else 0.0
    # Qg initial: use table QG
    for i, b in enumerate(barras):
        x0[idx['Qg'][0] + i] = barras_data[b]['QG']
    # Tap initial: use value from linhas list (we assume linha indices start at 1..nL)
    for ln in linhas:
        pos = idx['Tap'][0] + (ln['Linha'] - 1)
        x0[pos] = ln.get('Tap', 10.0)
    # s and r initial (positive)
    for i in range(nh):
        x0[idx['s'][0] + i] = 1.0
        x0[idx['r'][0] + i] = 1.0
    # lamb and pi initial zeros
    # already zero from initialization
    return x0

x0 = build_initial_x()

# --------------------------
# Bounds for variables
# --------------------------
bounds = []
# V bounds
for b in barras:
    bounds.append((Vmin_global, Vmax_global))
# Theta bounds: -pi..pi
for _ in barras:
    bounds.append((-math.pi, math.pi))
# Qg bounds
for b in barras:
    bounds.append((barras_data[b]['Qmin_G'], barras_data[b]['Qmax_G']))
# Tap bounds
for ln in linhas:
    bounds.append((Tap_min_global, Tap_max_global))  # fixed 10..10
# s bounds (positive)
for _ in range(nh):
    bounds.append((min_slack, None))
# r bounds (positive)
for _ in range(nh):
    bounds.append((min_slack, None))
# lambda bounds (unbounded)
for _ in range(ng):
    bounds.append((None, None))
# pi bounds (unbounded)
for _ in range(nh):
    bounds.append((None, None))

# --------------------------
# Otimização externa (loop de barreira/penalidade)
# --------------------------
mu_val = mu_init
rho_val = rho_init
x_current = x0.copy()

print("Iniciando otimização LPIP (numérica). Variáveis totais:", n_vars)
for outer in range(outer_iters):
    print(f"\n--- Iteração externa {outer+1}/{outer_iters} (mu={mu_val:.3e}, rho={rho_val:.3e}) ---")
    # wrapper objective for minimize (only needs x)
    fun = lambda x: objective_for_optimizer(x, mu_val, rho_val)

    # Run optimizer (SLSQP)
    res = minimize(
        fun=fun,
        x0=x_current,
        method='SLSQP',
        bounds=bounds,
        options={'maxiter': 100, 'ftol': tol_inner, 'disp': True}
    )

    if not res.success:
        print("Aviso: otimizador interno não convergiu:", res.message)
    else:
        print("Otimizador interno convergiu.")

    x_current = res.x.copy()

    # checar violação de restrições
    gx = g_vector(x_current)
    hx = h_vector(x_current)
    max_eq_violation = np.max(np.abs(gx)) if gx.size > 0 else 0.0
    max_ineq_violation = np.max(np.maximum(hx, 0.0)) if hx.size > 0 else 0.0
    print(f"Máx violação igualdade |g| = {max_eq_violation:.6e}, máx violação desigualdade h+ = {max_ineq_violation:.6e}")
    print("Valor objetivo f(x) (perdas):", f_loss(x_current))

    # Critério de parada simples
    if max_eq_violation < 1e-6 and max_ineq_violation < 1e-6 and mu_val < 1e-8:
        print("Convergiu externamente (critério de tolerância atendido).")
        break

    # Atualizar parâmetros
    mu_val *= mu_reduction
    # opcional: aumentar rho para reforçar objetivo, aqui mantemos igual
    # rho_val *= 1.0

# --------------------------
# Exibir resultados principais
# --------------------------
V_opt = get_slice(x_current, 'V')
Teta_opt = get_slice(x_current, 'Teta')
Qg_opt = get_slice(x_current, 'Qg')
Tap_opt = get_slice(x_current, 'Tap')

print("\n--- RESULTADOS FINAIS ---")
for i, b in enumerate(barras):
    print(f"Barra {b}: V = {V_opt[i]:.6f}, Teta = {Teta_opt[i]:.6f}, Qg = {Qg_opt[i]:.6f}")
print("\nTaps (por linha):")
for ln in linhas:
    pos = ln['Linha'] - 1
    print(f"Linha {ln['Linha']} ({ln['orig']}->{ln['dest']}): Tap = {Tap_opt[pos]:.6f}")

print("\nFunção objetivo final (perdas):", f_loss(x_current))
print("Máx violação igualdade final:", np.max(np.abs(g_vector(x_current))))
print("Máx violação desigualdade final:", np.max(np.maximum(h_vector(x_current), 0.0)))