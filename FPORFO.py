import pandas as pd
import numpy as np
import math
import mysql.connector

conexao = mysql.connector.connect( # conexão com banco de dados
  host='localhost',
  user='root',
  password='',
  database='sistema_eletrico'
)

cursor = conexao.cursor() # Criando cursor para fazer as solicitações ao banco

linha = input("Qual a linha que será calculada?\n") # Linha para a consulta no banco para poder realizar a função objetiva

cursor.execute("SELECT R FROM dadoslinha WHERE Linha = %s", (linha,)) # Solicitação do R em dados linha para calcular o Gkm
Rkm = cursor.fetchone() # Fetchone para armazenar apenas um dado
cursor.execute("SELECT X FROM dadoslinha WHERE Linha = %s", (linha,))
Xkm = cursor.fetchone()

if Rkm and Xkm: # Verifica se nenhum dos dois é nulo
    Rkm = Rkm[0]  # Pega os valores
    Xkm = Xkm[0]  

    print(f"Rkm: {Rkm}\nXkm: {Xkm}")

    Gkm = Rkm / (Rkm**2 + Xkm**2) # Calcula o Gkm
    Bkm = Xkm / (Rkm**2 + Xkm**2) # Calcula o Bkm para caso precise usar

    print(f"\nGkm: {Gkm}")
    print(f"\nBkm: {Bkm}")
    
else:
    print("Erro ao calcular o Gkm, um dos dados não existe")


cursor.execute("SELECT Tap FROM dadoslinha WHERE Linha = %s", (linha,))
Tkm = cursor.fetchone() # Tap da linha selecionada

Tkm = Tkm[0] # Puxar o valor

if Tkm == 0:
    Tkm = 1 # Caso seja 0, TAP será usado como 1
    

cursor.execute("SELECT V FROM dadosbarra WHERE Barra = (SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s)", (linha,))
Vk = cursor.fetchone() # V da Barra origem
cursor.execute("SELECT V FROM dadosbarra WHERE Barra = (SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s)", (linha,))
Vm = cursor.fetchone() # V da Barra Destino

Vk = Vk[0] # Puxar o Valor
Vm = Vm[0] # Puxar o Valor

cursor.execute("SELECT Theta FROM dadosbarra WHERE Barra = (SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s)", (linha,))
ThetaK = cursor.fetchone() # Ângulo da barra origem
cursor.execute("SELECT Theta FROM dadosbarra WHERE Barra = (SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s)", (linha,))
ThetaM = cursor.fetchone() # Ângulo da barra destino

ThetaK = ThetaK[0]
ThetaM = ThetaM[0]

ThetaKM = ThetaK - ThetaM # Diferença entre ambos os Ângulos

ThetaKM_rad = math.radians(ThetaKM) # Transformar em Radiano para fazer o Cosseno
cosThetaKM = math.cos(ThetaKM_rad) # Cosseno do ThetaKM (DIferença entre ambos os Ângulos)

print(f"Tkm: {Tkm}\nVk: {Vk}\nVm: {Vm}\nThetaK: {ThetaK}\nThetaM: {ThetaM}\nThetaKM: {ThetaKM}\nThetaKM_rad: {ThetaKM_rad}\nCosThetaKM: {cosThetaKM}") # Exibir todas variáveis para o cálculo

FuncaoObjetivo = Gkm*((1/(Tkm**2))*pow(Vk,2)+pow(Vm,2) - 2*(1/Tkm)*Vk*Vm*cosThetaKM) # Função Objetivo2


cursor.execute("SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s", (linha,))
BarraOrigem = cursor.fetchone() # Pegar a barra origem que estamos trabalhando
cursor.execute("SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s", (linha,))
BarraDestino = cursor.fetchone() # Pegar a Barra Destino que estamos trabalhando

BarraOrigem = BarraOrigem[0]
BarraDestino = BarraDestino[0]


print(f"\nFunção Objetiva da linha {linha} (barra {BarraOrigem} -> barra {BarraDestino}): {FuncaoObjetivo}") # Exibir a Função Objetivo


print(f"\nFunção Objetiva Escrita ({linha} (barra {BarraOrigem} -> barra {BarraDestino})): {Gkm}*((1/{Tkm}^2))*V{BarraOrigem}^2+V{BarraDestino}^2 - 2*(1/{Tkm})*V{BarraOrigem}*V{BarraDestino}*{cosThetaKM}") # Exibir a Função Objetivo1
