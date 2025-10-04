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

# =========================================================================
# FUNÇÃO OBJETIVO
# =========================================================================
for x in range(1, quantDeRegistros + 1):
    cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (x,))
    resultado = cursor.fetchone()

    if resultado:
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = resultado
        
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        Vk = result[0]
        
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        Vm = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        ThetaK = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        ThetaM = result[0]
        
        if Tkm == 0: Tkm = 1

        # Cálculo da Função Objetivo usando math (mais preciso que a subs, mas mantendo a subs para consistência da FO)
        resultado_termo = Gkm * (
            (1/Tkm**2) * Vk**2 + Vm**2 - 2 * (1/Tkm) * Vk * Vm * math.cos(ThetaK - ThetaM)
        )

        resultadoF += resultado_termo

        cosThetaKM = math.cos(ThetaK - ThetaM)
        termo = f"{Gkm}*((1/{Tkm}^2)*{Vk}^2 + {Vm}^2 - 2*(1/{Tkm})*{Vk}*{Vm}*{cosThetaKM}) \n"
        
        if x < quantDeRegistros:
            termo += " + "
        FuncaoObjetivoEscrita += termo

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

# Construção das restrições de potência ativa (Pkm)
funcRest1Escrita = ""

# Expressões para o fluxo de potência ativa (Pkm)
# CORREÇÃO: Pk para barra de origem (k)
def Pkm_Origem(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM):
    return (Gkm * (1/Tkm**2) * Vk**2) - ((1/Tkm) * Vk * Vm * (Gkm * math.cos(ThetaK - ThetaM) + Bkm * math.sin(ThetaK - ThetaM)))

# CORREÇÃO: Pmk para barra de destino (k)
def Pkm_Destino(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM):
    # O Vk aqui é a tensão da barra de destino (i) e Vm é a tensão da origem da linha
    return (Gkm * Vk**2) - ((1/Tkm) * Vk * Vm * (Gkm * math.cos(ThetaK - ThetaM) + Bkm * math.sin(ThetaK - ThetaM)))


for i in range(1, quantDeBarras + 1):
    if i == ignora:
        continue

    # Busca dados da barra i
    cursor.execute("SELECT Tipo, PG, QG, Qmin_G, Qmax_G, PC, QC, Bsh FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Pg, Qg, Qmin_G, Qmax_G, Pc, Qc, Bsh_barra = cursor.fetchone()

    # Variável acumuladora para o fluxo de potência ativa injetada na barra i
    potencia_fluxo_i = 0 
    termos = []
    
    # === Fluxo de potência Pkm, onde i é a Barra de ORIGEM (k) ===
    for linha in linhasPorBarra_Origem[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = cursor.fetchone()
        if Tkm == 0: Tkm = 1

        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        Vk = result[0]
        
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        Vm = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        ThetaK = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        ThetaM = result[0]

        cosThetaKM = math.cos(ThetaK - ThetaM)
        sinThetaKM = math.sin(ThetaK - ThetaM)

        # CÁLCULO DIRETO COM FLOAT (CORREÇÃO DE PRECISÃO)
        potencia_fluxo_linha = Pkm_Origem(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM)
        potencia_fluxo_i += potencia_fluxo_linha 
        
        termo = f"(({Gkm}*(1/{Tkm}**2))*({Vk}**2) - ((1/{Tkm})*{Vk})*{Vm}*((({Gkm})*({cosThetaKM})) + ({Bkm})*({sinThetaKM}))) \n"
        termos.append(termo)


    # === Fluxo de potência Pkm, onde i é a Barra de DESTINO (k) - Pmk ===
    for linha in linhasPorBarra_Destino[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = cursor.fetchone()
        
        if Tkm == 0: Tkm = 1
            
        # Para o fluxo M->K, V_k é a tensão da barra i (destino) e V_m é a da origem (BarraOrigem)
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        Vm_origem = result[0] 
        
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        Vk_destino = result[0] # V da barra i
        
        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        ThetaM_origem = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        ThetaK_destino = result[0] # Teta da barra i

        cosThetaKM = math.cos(ThetaK_destino - ThetaM_origem)
        sinThetaKM = math.sin(ThetaK_destino - ThetaM_origem)

        # CÁLCULO DIRETO COM FLOAT (CORREÇÃO DE PRECISÃO)
        potencia_fluxo_linha = Pkm_Destino(Gkm, Bkm, Tkm, Vk_destino, Vm_origem, ThetaK_destino, ThetaM_origem)
        potencia_fluxo_i += potencia_fluxo_linha
        
        # A string de visualização mantém a sua sintaxe original (que é o que resulta no -41.5121 na calculadora)
        termo = f"(({Gkm}*({Vk_destino}**2)) - ((1/{Tkm})*{Vk_destino})*{Vm_origem}*((({Gkm})*({cosThetaKM})) + ({Bkm})*({sinThetaKM}))) \n"
        termos.append(termo)
    
    # O resultado final da restrição é a soma dos fluxos menos a potência líquida injetada (PG - PC)
    resultadoRest1F = potencia_fluxo_i - Pg + Pc
    
    funcRest1 = " + ".join(termos)
    funcRest1Escrita += f"\nRestrição da Barra {i}:\n{funcRest1}  - {Pg} + {Pc} = {resultadoRest1F}\n ----------------------- X --------------------- \n"


# Construção das restrições de potência reativa (Qkm)
funcRest2Escrita = ""

# Expressões para o fluxo de potência reativa (Qkm)
def Qkm_Origem(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM, Bsh_linha):
    return ((-1 * (Bkm * (1/Tkm**2)) + Bsh_linha) * Vk**2) + ((1/Tkm) * Vk * Vm * (Bkm * math.cos(ThetaK - ThetaM) - Gkm * math.sin(ThetaK - ThetaM)))

def Qkm_Destino(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM, Bsh_linha):
    # Vk é a tensão da barra i (destino)
    return ((-1 * (Bkm + Bsh_linha)) * Vk**2) + ((1/Tkm) * Vk * Vm * (Bkm * math.cos(ThetaK - ThetaM) - Gkm * math.sin(ThetaK - ThetaM)))


for i in range(1, quantDeBarras + 1):
    if i == ignora:
        continue

    cursor.execute("SELECT Tipo, PG, QG, Qmin_G, Qmax_G, PC, QC, Bsh FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Pg, Qg, Qmin_G, Qmax_G, Pc, Qc, Bshk = cursor.fetchone()
    
    if Tipo != 0:
        continue
    
    potencia_refluxo_i = 0
    termos = []

    # === Fluxo de potência Qkm, onde i é a Barra de ORIGEM (k) ===
    for linha in linhasPorBarra_Origem[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino, Bsh FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino, Bsh_linha = cursor.fetchone()

        if Tkm == 0: Tkm = 1

        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        Vk = result[0]
        
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        Vm = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        ThetaK = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        ThetaM = result[0]

        cosThetaKM = math.cos(ThetaK - ThetaM)
        sinThetaKM = math.sin(ThetaK - ThetaM)

        # CÁLCULO DIRETO COM FLOAT (CORREÇÃO DE PRECISÃO)
        potencia_refluxo_linha = Qkm_Origem(Gkm, Bkm, Tkm, Vk, Vm, ThetaK, ThetaM, Bsh_linha)
        potencia_refluxo_i += potencia_refluxo_linha 

        termo = f"((-1*(({Bkm}*(1/{Tkm}**2))) + {Bsh_linha}) * ({Vk}**2) + ((1/{Tkm})*{Vk})*{Vm}*(({Bkm}*{cosThetaKM}) - ({Gkm}*{sinThetaKM}))) \n"
        termos.append(termo)
    
    # === Fluxo de potência Qmk, onde i é a Barra de DESTINO (k) ===
    for linha in linhasPorBarra_Destino[i]:

        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino, Bsh FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino, Bsh_linha = cursor.fetchone()

        if Tkm == 0: Tkm = 1

        # V_m é a tensão da barra de origem (BarraOrigem) e V_k é a da destino (BarraDestino = i)
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        Vm_origem = result[0] 
        
        cursor.execute("SELECT V FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        Vk_destino = result[0] 

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraOrigem,))
        result = cursor.fetchone()
        ThetaM_origem = result[0]

        cursor.execute("SELECT Teta FROM dadosbarra WHERE barra = %s", (BarraDestino,))
        result = cursor.fetchone()
        ThetaK_destino = result[0]

        cosThetaKM = math.cos(ThetaK_destino - ThetaM_origem)
        sinThetaKM = math.sin(ThetaK_destino - ThetaM_origem)

        # CÁLCULO DIRETO COM FLOAT (CORREÇÃO DE PRECISÃO)
        potencia_refluxo_linha = Qkm_Destino(Gkm, Bkm, Tkm, Vk_destino, Vm_origem, ThetaK_destino, ThetaM_origem, Bsh_linha)
        potencia_refluxo_i += potencia_refluxo_linha 

        termo = f"((-1*({Bkm} + {Bsh_linha})) * ({Vk_destino}**2) + ((1/{Tkm})*{Vk_destino})*{Vm_origem}*(({Bkm}*{cosThetaKM}) - ({Gkm}*{sinThetaKM}))) \n"
        termos.append(termo)

    # QSHk é a potência reativa do shunt da barra: Bshk * V_i^2
    QSHk = Bshk * (Vk_destino**2 if 'Vk_destino' in locals() else 0) 
    
    resultadoRest2F = potencia_refluxo_i - Qg + Qc - QSHk

    restricao2 = " + ".join(termos)
    funcRest2Escrita += f"\nRestrição da Barra {i}:\n{restricao2}  - {Qg} + {Qc} - {QSHk} = {resultadoRest2F}\n ----------------------- X --------------------- \n"


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