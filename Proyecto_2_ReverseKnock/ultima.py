import cobra

from gurobipy import GRB
from cobra import Reaction

# carga limpia del modelo(Usando el JSON) (si no ponemos esto se llena como de mensajes y se ve feo jaja)
model = cobra.io.load_json_model("macrofago_limpio.json")
model.solver.problem.setParam('NonConvex', 2)

# creamos atpm(Siempre necesaria porque el modelo base no traía esto)
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

# definimos las primeras cosas importantes:
biomasa_rxn = model.reactions.ATPM
target = model.reactions.EX_prostgd2_LPAREN_e_RPAREN_ 

#definimos constantes y una ayuda para mas adelante:

prob = model.problem #para ahorrarnos escribir model.problem en el futuro

M = 1000 #Para el big M

K = 5 #numero maximo de knockouts que vamos a permitir

f = 0.1  #es lo minimo que queremos que crezca la célula

restriccion_supervivencia = prob.Constraint(biomasa_rxn.flux_expression, lb=f,name="MinimaBiomasa")

model.add_cons_vars([restriccion_supervivencia])
#Primero vamos a restringir las reacciones que podemos borrar y las que no (para no interferir con la celula)
#porque para el programa es mas facil matar la celula asi que vamos a impedirle hacer eso:

validos = []

for r in model.reactions:
    #Intercambios, Biomasa, Target y Transportes (estas restricciones las saque del libro):
    if "EX_" not in r.id and r != biomasa_rxn and r != target:
        validos.append(r)

    #if EX_ not in r.id and r != biomasa_rxn and r != target:


#Definimos nuestra variable binaria (representa si la reaccion esta o no esta activa solo para los validos que definicmos hace poco)
y_vars={}
for r in validos:
    y = prob.Variable('y_' + r.id,type='binary')
    y_vars[r.id] = y

#ahora lo agregamos al modelo como variable:
model.add_cons_vars(list(y_vars.values()))

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
    mu_lb = dual_mu_lb[r.id]
    c_j=0

    #esto es lo que me va a ahorrar esto dos veces:
    if r==biomasa_rxn:
        c_j=1

    ecu=prob.Constraint(suma+mu_ub-mu_lb-c_j,lb=0,ub=0,name="EcuDual"+r.id)
    
    restriccion.append(ecu)
    
#Luego de definir la lista, la agregamos como restriccion:
model.add_cons_vars(restriccion)


restricciones_810_811 = []
for r in validos:
    y = y_vars[r.id]
    mu_u = dual_mu_ub[r.id]
    mu_l = dual_mu_lb[r.id]
    
    # Ecuación correcta: mu <= M * (1 - y)
    cons_810 = prob.Constraint(mu_u - M * (1 - y), ub=0, name='Dual_810_' + r.id)
    cons_811 = prob.Constraint(mu_l - M * (1 - y), ub=0, name='Dual_811_' + r.id)
    restricciones_810_811.extend([cons_810, cons_811])

model.add_cons_vars(restricciones_810_811)

#Esta es la parte final de definir cosas!
#Vamos a definir la restriccion de la dualidad fuerte

resultados = []

for r in model.reactions:
    # 1. Definir los límites (igual que antes)
    if r.upper_bound >= 1000:
        val_up = M
    else:
        val_up = r.upper_bound
    
    if r.lower_bound <= -1000:
        val_lb = -M
    else:
        val_lb = r.lower_bound

    mu_u = dual_mu_ub[r.id]
    mu_l = dual_mu_lb[r.id]

    # 2. Verificar si la reacción es "cortable"
    if r.id in y_vars:
        # CASO A: Es cortable (está en validos). Multiplicamos por la variable binaria 'y'.
        y = y_vars[r.id]
        # Aquí está la magia: conectamos el dual con el corte
        termino = (mu_u * val_up * y) - (mu_l * val_lb * y)
    else:
        # CASO B: No es cortable (siempre activa). Usamos la fórmula original.
        termino = (mu_u * val_up) - (mu_l * val_lb)

    resultados.append(termino)


suma_total = sum(resultados)
const_dual = prob.Constraint(biomasa_rxn.flux_expression - suma_total, lb=0, ub=0, name="DualidadFuerte")
model.add_cons_vars([const_dual])

# --- EJECUTAR ---

# Objetivo: minimizar el objetivo:
model.objective = prob.Objective(target.flux_expression, direction='min')

solution = model.optimize()
print(f"Estado: {solution.status}")
print(f" Flujo : {solution.objective_value}")
print(f" Flujo Biomasa: {biomasa_rxn.flux}")

cortes=[]
# Ver qué cortó
print("CORTES:")
for r_id, y in y_vars.items():
    if y.primal < 0.1: print(f"❌ {r_id}")
    cortes.append(r_id)