# Fluxo de Potência Ótimo Reativo (FPOR)  
### Optimal Reactive Power Flow (ORPF)

Este projeto implementa a formulação do Fluxo de Potência Ótimo Reativo (FPOR) em Python,  
utilizando MySQL como banco de dados para armazenar as barras e linhas do sistema elétrico.  
As restrições e a função objetivo são geradas automaticamente com base no trabalho do **Prof. Rafael Ramos de Souza (UNESP – Faculdade de Engenharia)**.  

---

This project implements the formulation of the Optimal Reactive Power Flow (ORPF) in Python,  
using MySQL as the database to store the buses and lines of the electric system.  
It automatically generates the objective function and constraints based on the work of **Prof. Rafael Ramos de Souza (UNESP – School of Engineering)**.  

---


---

## Tecnologias usadas / Technologies
- Python **3.13**
- MySQL
- Sympy
- Pandas
- Numpy  

---

---

## 🚀 Como executar / How to run

### Em Português
1. Configure o banco de dados **`sistema_eletrico`** com as tabelas `dadosbarra` e `dadoslinha`.  
2. Ajuste as credenciais no código Python.  
3. Execute o projeto:  
   ```bash
   python main.py

### In English
1. Set up the **`sistema_eletrico`** database with the tables `dadosbarra` and `dadoslinha`.
2. Configure the database credentials in the Python file.
3. Run the project:
    ```bash
     python main.py

## References
Rafael Ramos de Souza – UNESP, School of Engineering
Jéssica Antonio Delgado - UNESP, School of Engineering

## Obs...
- Este projeto é de **caráter acadêmico**.
- As fórmulas implementadas podem ser estendidas para outros métodos de otimização.

---
  
- This project is for **academic purposes**.  
- The implemented formulas can be extended to other optimization methods.
