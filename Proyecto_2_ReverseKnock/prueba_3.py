import cobra
import os
import sys
import contextlib
from cobra import Reaction

from gurobipy import GRB

# carga limpia del modelo(Usando el JSON) (si no ponemos esto se llena como de mensajes y se ve feo jaja)
model = cobra.io.load_json_model("macrofago_limpio.json")
model.solver="gurobi"
model.solver.problem.setParam('NonConvex', 2)

# creamos atpm (Siempre necesaria porque el modelo base no traía esto)
if "ATPM" not in model.reactions:
    atpm = Reaction('ATPM')
    atpm.name = 'ATP Maintenance'
    atpm.lower_bound = 0
    atpm.upper_bound = 1000
    atpm.add_metabolites({
        model.metabolites.atp_c: -1,
        model.metabolites.h2o_c: -1,
        model.metabolites.adp_c: 1,
        model.metabolites.pi_c:  1,
        model.metabolites.h_c:   1
    })
    model.add_reactions([atpm])
######################################################################################################
# definimos las primeras cosas importantes:
biomasa_rxn = model.reactions.ATPM

#ahora hay 2 targets asi que hice esto:
target_ldh=model.reactions.LDH_L
target_EX=model.reactions.EX_lac_L_LPAREN_e_RPAREN_
#pensaba volverlo una lista y hacer lo mismo pero para que se pueda hacer de forma mas general

#y este es el oficial:
target = target_ldh.flux_expression + target_EX.flux_expression


######################################################################################################
#definimos constantes y una ayuda para mas adelante:
prob = model.problem #para ahorrarnos escribir model.problem en el futuro

M = 3000 #Para el big M

K = 5 #numero maximo de knockouts que vamos a permitir
#Primero vamos a restringir las reacciones que podemos borrar y las que no (para no interferir con la celula)
#porque para el programa es mas facil matar la celula asi que vamos a impedirle hacer eso:


######################################################################################################
#estas son las reacciones que permitiremos cortar
validos = []
for r in model.reactions:
    #vamos a sacar la posibilidad de eliminar la biomasa
      if r != biomasa_rxn and r != target_ldh and r!=target_EX and "EX_" not in r.id:
        validos.append(r)
#######################################################################################################


#Definimos nuestra variable binaria (representa si la reaccion esta o no esta activa solo para los validos que definicmos hace poco)
y_vars={}
for r in validos:
    y = prob.Variable('y_' + r.id,type='binary')
    y_vars[r.id] = y
#ahora lo agregamos al modelo como variable:
model.add_cons_vars(list(y_vars.values()))

######################################################################################################

#Ahora vamos a empezar definiendo los Upper Bound (up) y los Lower Bound (lb) usando la tecnica del bigM:
cons_primal = [] #Estas van a ser los limites superiores e inferiores de cada reaccion en "validos"

terminos = [] #lista de los (1-y_j) para que la suma sea menor o igual a el valor de K de antes.

for r in validos:
    y = y_vars[r.id]
    
    upper_restriccion = prob.Constraint(r.flux_expression - M*y,ub=0,name='Big_M_Up'+r.id)
    cons_primal.append(upper_restriccion)
    
    if r.lower_bound <0: #esto pasa solo si la reaccion es reversible
        lower_restriccion = prob.Constraint(-M*y-r.flux_expression,ub=0,name='Big_M_low'+r.id)
        cons_primal.append(lower_restriccion)
    terminos.append(1 - y)

# agregamos como restriccion a los limites de la lista que definimos recien:
model.add_cons_vars(cons_primal)

######################################################################################################

#Además agregamos que la suma de los (1-y_j) sea menor o igual a K que es el numero de knockouts que estamos permitiendo
model.add_cons_vars(prob.Constraint(sum(terminos),ub=K,name='Max_Knockouts'))
#Ahora toca escribir las restricciones duales lambda y mu como diccionarios:
dual_lambda={}
dual_mu_ub={}
dual_mu_lb={}

for m in model.metabolites:
    nombre='lambda'+m.id 
    lam = prob.Variable(nombre,lb=-10000,ub=10000)
    dual_lambda[m.id]=lam

for r in model.reactions:
    nombre_ub='mu_ub'+r.id
    mu_u=prob.Variable(nombre_ub,lb=0,ub=10000)
    nombre_lb='mu_lb'+r.id
    mu_l=prob.Variable(nombre_lb,lb=0,ub=10000)
    dual_mu_lb[r.id]=mu_l
    dual_mu_ub[r.id]=mu_u

model.add_cons_vars(list(dual_lambda.values()))
model.add_cons_vars(list(dual_mu_lb.values()))
model.add_cons_vars(list(dual_mu_ub.values()))

######################################################################################################

#Ahora voy a escribir la suma rara de la matriz estequimetrica =1 con un if para no hacer esto dos veces

restriccion=[] 

for r in model.reactions: 
    terminos_suma=[]

    for m in r.metabolites:
        coefi=r.metabolites[m]
        lam=dual_lambda[m.id]
        terminos_suma.append(coefi*lam)

    suma=sum(terminos_suma)


    mu_ub = dual_mu_ub[r.id]
    mu_l = dual_mu_lb[r.id]
    
    c_j = 0
    
    if r == biomasa_rxn:
        c_j = 1
    elif r == target_ldg or r==target_EX:
        c_j = 0.0001 
        
    # El resto de la ecuación sigue igual
    ecu = prob.Constraint(suma + mu_ub - mu_l - c_j, lb=0.01, ub=0.01, name="EcuDual"+r.id)
    
    restriccion.append(ecu)

model.add_cons_vars(restriccion)
######################################################################################################
restricciones_nuevas = []
for r in validos:
    y = y_vars[r.id]
    mu_u = dual_mu_ub[r.id]
    mu_l = dual_mu_lb[r.id]
    
    # Ecuación correcta: mu <= M * (1 - y)
    cons_810 = prob.Constraint(mu_u - M * (1 - y), ub=0, name='Dual_810_' + r.id)
    cons_811 = prob.Constraint(mu_l - M * (1 - y), ub=0, name='Dual_811_' + r.id)
    restricciones_nuevas.extend([cons_810, cons_811])

model.add_cons_vars(restricciones_nuevas)

######################################################################################################


#Esta es la parte final de definir cosas!
#Vamos a definir la restriccion de la dualidad fuerte

# --- duali fuerte lineal ahora si!!!!

M_dual = 10000  # Debe ser igual o mayor al upper bound de tus variables duales

resultados_dualidad = []
restricciones_z = []

for r in model.reactions:
    # Definimos los límites!!
    if r.upper_bound >= 1000:
        val_up = M
    else:
        val_up = r.upper_bound
    
    if r.lower_bound <= -1000:
        val_lb = -M
    else:
        val_lb = r.lower_bound

    #Obtener las variables duales
    mu_u = dual_mu_ub[r.id]
    mu_l = dual_mu_lb[r.id]

    # 3. Verificar si la reacción tiene variable binaria 'y'
    if r.id in y_vars:
        y = y_vars[r.id]

    #creamos variables auxiliares (una para upper y otra para lower):
        z_u=prob.Variable(f"z_up_{r.id}",lb=0,ub=M_dual)
        z_l=prob.Variable(f"z_l_{r.id}", lb=0,ub=M_dual)
    
    #ahora tenemos 3 restricciones para cada una:
    #Este para upper
        restricciones_z.append(prob.Constraint(z_u - M_dual*y, ub=0))
        restricciones_z.append(prob.Constraints(z_u-mu_u, ub=0))
        restricciones_z.append(prob.Constraints(mu_u - z_u - M_dual*(1-y),ub=0))

    #este para lower:
        restricciones_z.append(prob.Constraint(z_l - M_dual*y, ub=0))
        restricciones_z.append(prob.Constraints(z_l-mu_l, ub=0))
        restricciones_z.append(prob.Constraints(mu_l- z_l - M_dual*(1-y),ub=0))

        termino=(z_u*val_up)-(z_l*val_lb)
    else:
        termino=(mu_u*val_up)-(mu_l*val_ub)
    resultados_dualidad.append(termino)

# Agregar todas las nuevas restricciones al modelo
model.add_cons_vars(restricciones_z)

# Suma total y restricción final
suma_total = sum(resultados_dualidad) 

coeficiente_inclinacion = 0.0001 
objetivo_primal_real = biomasa_rxn.flux_expression + (coeficiente_inclinacion)*(target_EX.flux_expression+target_ldh.flux_expression)

# Actualizamos la restricción de igualdad fuerte
const_dual = prob.Constraint(objetivo_primal_real- suma_total, lb=0, ub=0, name="DualidadFuerte")

model.add_cons_vars([const_dual])
######################################################################################################


f = 0.1  #Esto es cuanto queremos que la celula crezca como minimo (para que no se vaya a 0 le decimos que tiene que ser un f% de lo que es originalmente)

# 2. Agregamos la restricción: "Prohibido morir"
restriccion_supervivencia = prob.Constraint(biomasa_rxn.flux_expression, lb=f,name="MinimaBiomasa")

model.add_cons_vars([restriccion_supervivencia])


######################################################################################################

#ejecutamos todo y esperamos por lo mejor

model.objective = prob.Objective(target.flux_expression, direction='min')

solution = model.optimize()
print(f" Flujo del objetivo : {solution.objective_value}")
print(f" Flujo de la Biomasa: {biomasa_rxn.flux}")


#para que nos muestre que cosas cortó:

print("Genes cortados:")
for r_id, y in y_vars.items():
    if y.primal < 0.1: 
        print(f" {r_id}") 
