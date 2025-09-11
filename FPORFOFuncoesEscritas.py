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


# Definindo variáveis simbólicas fora do loop
GkmR, TkmR, VkR, VmR, ThetaKR, ThetaMR = sp.symbols("Gkm Tkm Vk Vm ThetaK ThetaM")

# Expressão simbólica da função objetivo
funcaoResol = "Gkm*((1/Tkm**2)*Vk**2 + Vm**2 - 2*(1/Tkm)*Vk*Vm*cos(ThetaK - ThetaM))"
expr = sp.sympify(funcaoResol, locals={"cos": sp.cos})


resultadoF = 0
FuncaoObjetivoEscrita = ""


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

        cosThetaKM = math.cos(ThetaK - ThetaM)

        resultado = expr.subs({
            GkmR: Gkm,
            TkmR: Tkm if Tkm != 0 else 1,
            VkR: Vk,
            VmR: Vm,
            ThetaKR: ThetaK,
            ThetaMR: ThetaM
        })

        resultadoF += resultado

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
    if origem == ignora or destino == ignora:
        continue
    linhasPorBarra_Origem[origem].append(linha)
    linhasPorBarra_Destino[destino].append(linha)

# Construção das restrições de potência ativa (Pkm)
funcRest1Escrita = ""
funcRest1 = ""
for i in range(1, quantDeBarras + 1):
    if i == ignora:
        continue

    cursor.execute("SELECT Tipo, PG, QG, Qmin_G, Qmax_G, PC, QC, Bsh FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Pg, Qg, Qmin_G, Qmax_G, Pc, Qc, Bsh = cursor.fetchone()

    termos = []

    for linha in linhasPorBarra_Origem[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = cursor.fetchone()
        if Tkm == 0:
            Tkm = 1

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

        termo = f"(({Gkm}*(1/{Tkm}**2))*({Vk}**2) - ((1/{Tkm})*{Vk})*{Vm}*((({Gkm})*({cosThetaKM})) + ({Bkm})*({sinThetaKM})) \n"
        termos.append(termo)

    for linha in linhasPorBarra_Destino[i]:
        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino = cursor.fetchone()

        if Tkm == 0:
            Tkm = 1


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

        termo = f"(({Gkm}*({Vk}**2)) - ((1/{Tkm})*{Vk})*{Vm}*((({Gkm})*({cosThetaKM})) + ({Bkm})*({sinThetaKM})) \n"
        termos.append(termo)
        
    funcRest1 = " + ".join(termos)
    funcRest1Escrita += f"\nRestrição da Barra {i}:\n{funcRest1}  - {Pg} + {Pc}\n ----------------------- X --------------------- \n"


funcRest2Escrita = ""

for i in range(1, quantDeBarras + 1):
    if i == ignora:
        continue

    cursor.execute("SELECT Tipo, PG, QG, Qmin_G, Qmax_G, PC, QC, Bsh FROM dadosbarra WHERE barra = %s", (i,))
    Tipo, Pg, Qg, Qmin_G, Qmax_G, Pc, Qc, Bshk = cursor.fetchone()
    if Tipo != 0:
        continue
    
    termos = []

    for linha in linhasPorBarra_Origem[i]:
        if Tipo != 0:
            continue

        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino, Bsh FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino, Bsh = cursor.fetchone()

        if Tkm == 0:
            Tkm = 1

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

        termo = f"(-1*(({Bkm}*(1/{Tkm}**2))) + {Bsh}) * ({Vk}**2) + ((1/{Tkm})*{Vk})*{Vm}*(({Bkm}*{cosThetaKM}) - ({Gkm}*{sinThetaKM})) \n"
        termos.append(termo)
    
    for linha in linhasPorBarra_Destino[i]:
        if Tipo != 0:
            continue

        cursor.execute("SELECT Gkm, Bkm, Tap, Barra_Origem, Barra_Destino, Bsh FROM dadoslinha WHERE Linha = %s", (linha,))
        Gkm, Bkm, Tkm, BarraOrigem, BarraDestino, Bsh = cursor.fetchone()

        if Tkm == 0:
            Tkm = 1

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

        termo = f"(-1*({Bkm} + {Bsh})) * ({Vk}**2) + ((1/{Tkm})*{Vk})*{Vm}*(({Bkm}*{cosThetaKM}) - ({Gkm}*{sinThetaKM})) \n"
        termos.append(termo)

    restricao2 = " + ".join(termos)
    funcRest2Escrita += f"\nRestrição da Barra {i}:\n{restricao2}  - {Qg} + {Qc} - 0\n ----------------------- X --------------------- \n"


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