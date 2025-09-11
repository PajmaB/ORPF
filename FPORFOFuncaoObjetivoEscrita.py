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

cursor.execute("SELECT * FROM dadoslinha")
registros = cursor.fetchall()

quantDeRegistros = 0
for i, registros in enumerate(registros, start=1):
   quantDeRegistros = i

FuncaoObjetivoEscrita = ""
x = 1 # Define a linha que está sendo interagido para escrever a funçao objetivo
while x >= 1:
 
 cursor.execute("SELECT R FROM dadoslinha WHERE Linha = %s", (x,)) # Solicitação do R em dados linha para calcular o Gkm
 Rkm = cursor.fetchone() # Fetchone para armazenar apenas um dado
 cursor.execute("SELECT X FROM dadoslinha WHERE Linha = %s", (x,))
 Xkm = cursor.fetchone()

 if Rkm and Xkm: # Verifica se nenhum dos dois é nulo
    Rkm = Rkm[0]  # Pega os valores
    Xkm = Xkm[0]  

    #print(f"Rkm: {Rkm}\nXkm: {Xkm}")

    Gkm = Rkm / (Rkm**2 + Xkm**2) # Calcula o Gkm
    Bkm = Xkm / (Rkm**2 + Xkm**2) # Calcula o Bkm para caso precise usar

    #print(f"\nGkm: {Gkm}")
    #print(f"\nBkm: {Bkm}")
    
 #else:
    
    #print("Erro ao calcular o Gkm, um dos dados não existe")

 cursor.execute("SELECT Tap FROM dadoslinha WHERE Linha = %s", (x,))
 Tkm = cursor.fetchone() # Tap da linha selecionada

 Tkm = Tkm[0] # Puxar o valor

 if Tkm == 0:
    Tkm = 1 # Caso seja 0, TAP será usado como 1
      

 cursor.execute("SELECT V FROM dadosbarra WHERE Barra = (SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s)", (x,))
 Vk = cursor.fetchone() # V da Barra origem
 cursor.execute("SELECT V FROM dadosbarra WHERE Barra = (SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s)", (x,))
 Vm = cursor.fetchone() # V da Barra Destino

 Vk = Vk[0] # Puxar o Valor
 Vm = Vm[0] # Puxar o Valor

 cursor.execute("SELECT Theta FROM dadosbarra WHERE Barra = (SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s)", (x,))
 ThetaK = cursor.fetchone() # Ângulo da barra origem
 cursor.execute("SELECT Theta FROM dadosbarra WHERE Barra = (SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s)", (x,))
 ThetaM = cursor.fetchone() # Ângulo da barra destino

 ThetaK = ThetaK[0]
 ThetaM = ThetaM[0]

 ThetaKM = ThetaK - ThetaM # Diferença entre ambos os Ângulos

 ThetaKM_rad = math.radians(ThetaKM) # Transformar em Radiano para fazer o Cosseno
 cosThetaKM = math.cos(ThetaKM_rad) # Cosseno do ThetaKM (Diferença entre ambos os Ângulos)
 sinThetaKM = math.sin(ThetaKM_rad)   # Seno do ThetaKM

 #print(f"Tkm: {Tkm}\nVk: {Vk}\nVm: {Vm}\nThetaK: {ThetaK}\nThetaM: {ThetaM}\nThetaKM: {ThetaKM}\nThetaKM_rad: {ThetaKM_rad}\nCosThetaKM: {cosThetaKM}") # Exibir todas variáveis para o cálculo

 FuncaoObjetivo = Gkm*((1/(Tkm**2))*pow(Vk,2)+pow(Vm,2) - 2*(1/Tkm)*Vk*Vm*cosThetaKM) # Função Objetivo

 cursor.execute("SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s", (x,))
 BarraOrigem = cursor.fetchone() # Pegar a barra origem que estamos trabalhando
 cursor.execute("SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s", (x,))
 BarraDestino = cursor.fetchone() # Pegar a Barra Destino que estamos trabalhando

 BarraOrigem = BarraOrigem[0]
 BarraDestino = BarraDestino[0]
 if x < quantDeRegistros:
    FuncLinha = f"{Gkm}*((1/{Tkm}^2))*V{BarraOrigem}^2+V{BarraDestino}^2 - 2*(1/{Tkm})*V{BarraOrigem}*V{BarraDestino}*{cosThetaKM} + "
    FuncaoObjetivoEscrita = FuncaoObjetivoEscrita + FuncLinha
 if x == quantDeRegistros:
    FuncLinha = f"{Gkm}*((1/{Tkm}^2))*V{BarraOrigem}^2+V{BarraDestino}^2 - 2*(1/{Tkm})*V{BarraOrigem}*V{BarraDestino}*{cosThetaKM}"
    FuncaoObjetivoEscrita = FuncaoObjetivoEscrita + FuncLinha
 x = x + 1
 if x == (quantDeRegistros + 1):
    x = 0
    
print(f"{FuncaoObjetivoEscrita}")

ignora = input("\nQual linha é a da barra Slack? \n") # Linha que será ignorada por ser a slack

x = 1

if(x == ignora): x = x + 1

while(x < quantDeRegistros): # Ira realizar os PKM
   
 cursor.execute("SELECT R FROM dadoslinha WHERE Linha = %s", (x,)) # Solicitação do R em dados linha para calcular o Gkm
 Rkm = cursor.fetchone() # Fetchone para armazenar apenas um dado
 cursor.execute("SELECT X FROM dadoslinha WHERE Linha = %s", (x,))
 Xkm = cursor.fetchone()

 if Rkm and Xkm: # Verifica se nenhum dos dois é nulo
    Rkm = Rkm[0]  # Pega os valores
    Xkm = Xkm[0]  

    #print(f"Rkm: {Rkm}\nXkm: {Xkm}")

    Gkm = Rkm / (Rkm**2 + Xkm**2) # Calcula o Gkm
    Bkm = Xkm / (Rkm**2 + Xkm**2) # Calcula o Bkm para caso precise usar

    #print(f"\nGkm: {Gkm}")
    #print(f"\nBkm: {Bkm}")
    
 #else:
    
    #print("Erro ao calcular o Gkm, um dos dados não existe")
    
 cursor.execute("SELECT Tap FROM dadoslinha WHERE Linha = %s", (x,))
 Tkm = cursor.fetchone() # Tap da linha selecionada

 Tkm = Tkm[0] # Puxar o valor

 if Tkm == 0:
    Tkm = 1 # Caso seja 0, TAP será usado como 1
      

 cursor.execute("SELECT V FROM dadosbarra WHERE Barra = (SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s)", (x,))
 Vk = cursor.fetchone() # V da Barra origem
 cursor.execute("SELECT V FROM dadosbarra WHERE Barra = (SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s)", (x,))
 Vm = cursor.fetchone() # V da Barra Destino

 Vk = Vk[0] # Puxar o Valor
 Vm = Vm[0] # Puxar o Valor

 cursor.execute("SELECT Theta FROM dadosbarra WHERE Barra = (SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s)", (x,))
 ThetaK = cursor.fetchone() # Ângulo da barra origem
 cursor.execute("SELECT Theta FROM dadosbarra WHERE Barra = (SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s)", (x,))
 ThetaM = cursor.fetchone() # Ângulo da barra destino

 ThetaK = ThetaK[0]
 ThetaM = ThetaM[0]

 ThetaKM = ThetaK - ThetaM # Diferença entre ambos os Ângulos

 ThetaKM_rad = math.radians(ThetaKM) # Transformar em Radiano para fazer o Cosseno
 cosThetaKM = math.cos(ThetaKM_rad) # Cosseno do ThetaKM (Diferença entre ambos os Ângulos)
 sinThetaKM = math.sin(ThetaKM_rad)   # Seno do ThetaKM

 cursor.execute("SELECT Barra_Origem FROM dadoslinha WHERE Linha = %s", (x,))
 BarraOrigem = cursor.fetchone() # Pegar a barra origem que estamos trabalhando
 cursor.execute("SELECT Barra_Destino FROM dadoslinha WHERE Linha = %s", (x,))
 BarraDestino = cursor.fetchone() # Pegar a Barra Destino que estamos trabalhando

 PkmInicial = Gkm*(1/pow(Tkm, 2))*pow(Vk, 2) - (1/Tkm)*Vk*Vm*(Gkm*cosThetaKM + Bkm*sinThetaKM)
 PkmFinal = Gkm*pow(Vk, 2) - ((1/Tkm)*Vk)*Vm*(Gkm*cosThetaKM + Bkm * sinThetaKM)
 Pkm = PkmInicial + PkmFinal
 
 print(f"Restrição 1 da linha {x} :)")
 x = x + 1