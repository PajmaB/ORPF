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