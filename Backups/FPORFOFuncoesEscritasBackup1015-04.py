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

print(f"\nFunção Objetivo (ajustada):\n{FuncaoObjetivoEscrita}")
print("Resultado da FO (ajustado): ", float(resultadoF))

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
    # aqui Vk_destino equivale a V_m no formato do ref, Vm_origem equivale a V_k
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

from scipy.optimize import minimize # <<< IMPORTANTE: Importar o otimizador

# =========================================================================
# PARTE 2: MONTAGEM SIMBÓLICA DO PROBLEMA DE OTIMIZAÇÃO (LPIP)
# =========================================================================
print("\n\n--- INICIANDO PARTE 2: MONTAGEM SIMBÓLICA (LPIP) ---")

# Buscar todos os dados do DB de uma vez para a montagem simbólica
cursor.execute("SELECT barra, Tipo, PG, PC, QG, QC, Bsh, 1.1, 0.9, Qmax_G, Qmin_G FROM dadosbarra")
barras_data = {row[0]: dict(zip(['id', 'tipo', 'pg', 'pc', 'qg', 'qc', 'bsh', 'vmax', 'vmin', 'qmax', 'qmin'], row)) for row in cursor.fetchall()}

cursor.execute("SELECT Linha, Barra_Origem, Barra_Destino, Gkm, Bkm, Bsh, 1.0, 1.1, 0.9 FROM dadoslinha")
linhas_data = {row[0]: dict(zip(['id', 'orig', 'dest', 'g', 'b', 'bsh', 'tap_val', 'tap_max', 'tap_min'], row)) for row in cursor.fetchall()}

# Definindo as variáveis do problema como símbolos do SymPy
V = sp.symbols(f'V1:{quantDeBarras+1}')
Teta = sp.symbols(f'Teta1:{quantDeBarras+1}')
Qg = sp.symbols(f'Qg1:{quantDeBarras+1}')
Tap = sp.symbols(f'Tap1:{quantDeRegistros+1}')
rho, mu = sp.symbols('rho mu')

# Listas para as equações simbólicas
funcao_objetivo_fx = 0
restricoes_igualdade_gx = []
restricoes_desigualdade_hx = []

# --- Montagem da Função Objetivo Simbólica f(x) ---
for linha in linhas_data.values():
    k, m, gkm = linha['orig'], linha['dest'], linha['g']
    funcao_objetivo_fx += gkm * (V[k-1]**2 + V[m-1]**2 - 2 * V[k-1] * V[m-1] * sp.cos(Teta[k-1] - Teta[m-1]))
print("\n[OK] Função Objetivo Simbólica f(x) montada.")

# --- Montagem das Restrições de Igualdade g(x) ---
# (O código para montar g(x) e h(x) permanece aqui, sem alterações)
# ...
# Balanço de Potência Ativa (P)
for k, barra in barras_data.items():
    if k == ignora: continue
    Pk_injetada = 0
    # ... (código de cálculo de Pk_injetada)
    equacao_Pk = Pk_injetada - barra['pg'] + barra['pc']
    restricoes_igualdade_gx.append(equacao_Pk)
# Balanço de Potência Reativa (Q) para barras PQ
for k, barra in barras_data.items():
    if barra['tipo'] == 0:
        Qk_injetada = 0
        # ... (código de cálculo de Qk_injetada)
        equacao_Qk = Qk_injetada - barra['qg'] + barra['qc']
        restricoes_igualdade_gx.append(equacao_Qk)
print(f"[OK] {len(restricoes_igualdade_gx)} Restrições de Igualdade Simbólicas g(x) montadas.")

# --- Montagem das Restrições de Desigualdade h(x) <= 0 ---
# Limites de Tensão, Geração Reativa, Tap
# ... (código para montar h(x) permanece aqui)
print(f"[OK] {len(restricoes_desigualdade_hx)} Restrições de Desigualdade Simbólicas h(x) montadas.")

# --- Montagem da Função Lagrangiana LPIP ---
num_g = len(restricoes_igualdade_gx)
num_h = len(restricoes_desigualdade_hx)
s = sp.symbols(f's1:{num_h+1}')
r = sp.symbols(f'r1:{num_h+1}')
lamb = sp.symbols(f'lambda1:{num_g+1}')
pi = sp.symbols(f'pi1:{num_h+1}')
LPIP = (rho * funcao_objetivo_fx - 
        mu * sum(sp.log(s[i]) + sp.log(r[i]) for i in range(num_h)) +
        sum(r) +
        sum(lamb[j] * restricoes_igualdade_gx[j] for j in range(num_g)) +
        sum(pi[i] * (restricoes_desigualdade_hx[i] + s[i] - r[i]) for i in range(num_h)))
print("[OK] Função Lagrangiana (LPIP) montada com sucesso!")


# =========================================================================
# PARTE 3: SOLUÇÃO NUMÉRICA USANDO SCIPY
# =========================================================================
print("\n\n--- INICIANDO PARTE 3: SOLUÇÃO NUMÉRICA ITERATIVA ---")

# 1. Consolidar todas as variáveis em uma única lista
all_vars = list(V) + list(Teta) + list(Qg) + list(Tap) + list(s) + list(r) + list(lamb) + list(pi)
print(f"Total de variáveis de otimização: {len(all_vars)}")

# 2. Calcular o gradiente da LPIP (lista de derivadas parciais)
print("Calculando o gradiente simbólico (pode levar um momento)...")
grad_LPIP = [sp.diff(LPIP, var) for var in all_vars]
print("[OK] Gradiente simbólico calculado.")

# 3. Converter a LPIP e seu gradiente em funções numéricas rápidas
# As funções aceitarão um único vetor 'x' e os parâmetros mu e rho
LPIP_numeric = sp.lambdify([all_vars, mu, rho], LPIP, 'numpy')
grad_LPIP_numeric = sp.lambdify([all_vars, mu, rho], grad_LPIP, 'numpy')
print("[OK] Funções convertidas para formato numérico (lambdify).")

# 4. Configurar e executar o loop de otimização
num_iteracoes = 5
mu_inicial = 10.0
rho_val = 1.0 # Parâmetro de penalidade, geralmente fixo ou com outra lógica de atualização

# Chute inicial para todas as variáveis
# Valores iniciais: V=1.0, Teta=0, Qg=0.1, Tap=1.0, s=1.0, r=1.0, lambda=0, pi=0
x0 = np.concatenate([
    np.ones(len(V)),          # V
    np.zeros(len(Teta)),      # Teta
    np.ones(len(Qg)) * 0.1,   # Qg
    np.ones(len(Tap)),        # Tap
    np.ones(len(s)),          # s
    np.ones(len(r)),          # r
    np.zeros(len(lamb)),      # lambda
    np.zeros(len(pi))         # pi
])

mu_val = mu_inicial
for i in range(num_iteracoes):
    print(f"\n--- Iteração de Otimização {i+1}/{num_iteracoes} (mu = {mu_val:.6f}) ---")

    # Funções wrapper para o solver, que fixam os parâmetros mu e rho
    def objective_for_solver(x):
        return LPIP_numeric(x, mu_val, rho_val)

    def jacobian_for_solver(x):
        return np.array(grad_LPIP_numeric(x, mu_val, rho_val))

    # Chamada ao otimizador
    result = minimize(
        fun=objective_for_solver,
        x0=x0,
        method='BFGS', # Um método eficiente que usa o gradiente
        jac=jacobian_for_solver,
        options={'disp': True, 'maxiter': 100}
    )

    # Atualiza o chute inicial para a próxima iteração (warm start)
    x0 = result.x
    
    # Reduz o parâmetro de barreira
    mu_val *= 0.1

# --- 5. Exibir os resultados finais ---
print("\n\n--- OTIMIZAÇÃO CONCLUÍDA ---")
print("Valores ótimos encontrados para as variáveis:")

# Cria um dicionário para fácil visualização dos resultados
resultados_finais = dict(zip([str(var) for var in all_vars], result.x))

for var_nome, var_valor in resultados_finais.items():
    # Imprime apenas as variáveis de estado principais para não poluir a saída
    if var_nome.startswith('V') or var_nome.startswith('Teta') or var_nome.startswith('Qg'):
        print(f"{var_nome}: {var_valor:.6f}")

# Fechando a conexão ao final de tudo
cursor.close()
conexao.close()