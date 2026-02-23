import cobra
import pandas as pd
from cobra import Reaction

#esto es como igual
model = cobra.io.load_json_model("macrofago_limpio.json")
model.solver = "gurobi"

# los parametros de antes 
M = 3000        
M_dual = 10000  
K = 5           
f = 0.1     

# hacemos esto antes

# ahora definimos nuestra solucion al problema como una funcion que vamos a usar mas abajo

def correr_optknock(genes_prohibidos, n):
    # En lugar de context manager, creamos una copia física nueva del modelo
    m = cobra.io.load_json_model("macrofago_limpio.json")
    prob = m.problem

    if "ATPM" not in m.reactions:
        atpm = Reaction('ATPM')
        atpm.name = 'ATP Maintenance'
        atpm.lower_bound = 0
        atpm.upper_bound = 1000
        atpm.add_metabolites({
            m.metabolites.atp_c: -1,
            m.metabolites.h2o_c: -1,
            m.metabolites.adp_c: 1,
            m.metabolites.pi_c:  1,
            m.metabolites.h_c:   1
        })
        m.add_reactions([atpm])
    
    # Identificamos las reacciones en la copia
    biomasa_rxn = m.reactions.ATPM
    target_ldh = m.reactions.LDH_L
    target_EX= m.reactions.EX_lac_L_LPAREN_e_RPAREN_
    target = target_ldh.flux_expression + target_EX.flux_expression
    
 ################################################################################################
    validos = []
    for r in m.reactions:
        # Filtro de reacciones que se pueden cortar
        if r != biomasa_rxn and r != target_ldh and r != target_EX and "EX_" not in r.id:
            if r.id not in genes_prohibidos:    #agregamos esto para no usar los genes que ya sacamos
                validos.append(r)

 #################################################################################################
    # Las BINARIAS:
    y_vars={}
    for r in validos:
        y = prob.Variable('y_' + r.id,type='binary')
        y_vars[r.id] = y
    #ahora lo agregamos al modelo como variable:
    m.add_cons_vars(list(y_vars.values()))

    target_export=0

 ################################################################################################
    # Restricciones Primales Big-M
    cons_primal = []
    terminos=[]
    for r in validos:
        y = y_vars[r.id]
        upper_restriccion = prob.Constraint(r.flux_expression - M*y,ub=0,name='Big_M_Up'+r.id)
        cons_primal.append(upper_restriccion)

        if r.lower_bound < 0:
            lower_restriccion = prob.Constraint(-M*y-r.flux_expression,ub=0,name='Big_M_low'+r.id)
            cons_primal.append(lower_restriccion)
        terminos.append(1 - y)    
    m.add_cons_vars(cons_primal)
    m.add_cons_vars(prob.Constraint(sum(terminos),ub=K,name='Max_Knockouts'))

 ################################################################################################
    # Variables Duales
    dual_lambda={}
    dual_mu_ub={}
    dual_mu_lb={}

    for met in m.metabolites:
        nombre='lambda'+met.id 
        lam = prob.Variable(nombre,lb=-10000,ub=10000)
        dual_lambda[met.id]=lam

    for r in m.reactions:
        nombre_ub='mu_ub'+r.id
        mu_u=prob.Variable(nombre_ub,lb=0,ub=10000)
        nombre_lb='mu_lb'+r.id
        mu_l=prob.Variable(nombre_lb,lb=0,ub=10000)
        dual_mu_lb[r.id]=mu_l
        dual_mu_ub[r.id]=mu_u

    m.add_cons_vars(list(dual_lambda.values()))
    m.add_cons_vars(list(dual_mu_lb.values()))
    m.add_cons_vars(list(dual_mu_ub.values()))

 #################################################################################################


    # Restricciones de Igualdad Dual
    duales = []
    for r in m.reactions: 
        terminos_suma=[]

        for met in r.metabolites:
            coefi=r.metabolites[met]
            lam=dual_lambda[met.id]
            terminos_suma.append(coefi*lam)

        suma=sum(terminos_suma)

        mu_ub = dual_mu_ub[r.id]
        mu_l = dual_mu_lb[r.id]
    
        c_j = 0
    
        if r == biomasa_rxn:
            c_j = 1
        elif r == target_ldh or r == target_EX:
            c_j = 0.0001 
        
    # El resto de la ecuación sigue igual
        ecu = prob.Constraint(suma + mu_ub - mu_l - c_j, lb=0, ub=0, name="EcuDual"+r.id)
    
        duales.append(ecu)

    m.add_cons_vars(duales)

 #################################################################################################

    restricciones_nuevas = []
    for r in validos:
        y = y_vars[r.id]
        mu_u = dual_mu_ub[r.id]
        mu_l = dual_mu_lb[r.id]
    
    # Ecuación correcta: mu <= M * (1 - y)
        cons_810 = prob.Constraint(mu_u - M * (1 - y), ub=0, name='Dual_810_' + r.id)
        cons_811 = prob.Constraint(mu_l - M * (1 - y), ub=0, name='Dual_811_' + r.id)
        restricciones_nuevas.extend([cons_810, cons_811])

    m.add_cons_vars(restricciones_nuevas)

 #################################################################################################


    # Dualidad Fuerte Linealizada!!
    terminos_duales = []
    restricciones_z = []
    for r in m.reactions:
    # 1. Definir los límites numéricos (igual que antes)
        if r.upper_bound >= M:
            val_up = M
        else:
            val_up = r.upper_bound
    
        if r.lower_bound <= -M:
            val_lb = -M
        else:
            val_lb = r.lower_bound

    # 2. Obtener las variables duales
        mu_u = dual_mu_ub[r.id]
        mu_l = dual_mu_lb[r.id]

        if r.id in y_vars:
            y = y_vars[r.id]
            z_u = prob.Variable(f"z_u_{r.id}_it{n}", lb=0, ub=M_dual)
            z_l = prob.Variable(f"z_l_{r.id}_it{n}", lb=0, ub=M_dual)
            
            restricciones_z.extend([
                prob.Constraint(z_u - M_dual * y, ub=0, name=f'lin_u1_{r.id}_it{n}'),
                prob.Constraint(z_u - mu_u, ub=0, name=f'lin_u2_{r.id}_it{n}'),
                prob.Constraint(mu_u - z_u - M_dual * (1 - y), ub=0, name=f'lin_u3_{r.id}_it{n}'),

                prob.Constraint(z_l - M_dual * y, ub=0, name=f'lin_l1_{r.id}_it{n}'),
                prob.Constraint(z_l - mu_l, ub=0, name=f'lin_l2_{r.id}_it{n}'),
                prob.Constraint(mu_l - z_l - M_dual * (1 - y), ub=0, name=f'lin_l3_{r.id}_it{n}')
            ])
            termino = (z_u * val_up) - (z_l * val_lb)
        else:
            termino = (mu_u * val_up) - (mu_l * val_lb)
        terminos_duales.append(termino)

    m.add_cons_vars(restricciones_z)
    
    coef_inc = 0.0001
    objetivo_primal_real = biomasa_rxn.flux_expression + coef_inc * (target_EX.flux_expression + target_ldh.flux_expression)
    suma_total=sum(terminos_duales)
 # Actualizamos la restricción de igualdad fuerte con un margen más generoso
    const_dual = prob.Constraint(objetivo_primal_real - suma_total, lb=-0.1, ub=0.1, name="DualidadFuerte")
    m.add_cons_vars([const_dual])
    # Supervivencia
    restriccion_supervivencia = prob.Constraint(biomasa_rxn.flux_expression, lb=f,name="MinimaBiomasa")
    m.add_cons_vars([restriccion_supervivencia])    
    # Objetivo: Minimizar Lactato

    m.objective = prob.Objective(target, direction='min')
    solution = m.optimize()
    
    #ahora con esto vemos el return de nuestra funcion:
    if solution.status == 'optimal':
        cortes = []
        for r_id, y in y_vars.items():
            if y.primal < 0.1: 
                cortes.append(r_id) # Guárdalo en la lista de cortados
        return cortes, solution.objective_value, biomasa_rxn.flux
    else:
        return [], None, None

#ahora esta es la pate de la "tierlist"
tier_list = []
acumulados = []

for i in range(5):
    print("Iteración "+str({i+1}))
    cortes, t_flux, b_flux = correr_optknock(acumulados, i + 1)
    
    if not cortes: 
        print("No se encontraron más soluciones óptimas.")
        break
        
    print(f"Cortes encontrados: {cortes}")
    tier_list.append({"Iteracion": i+1, "Lactato": t_flux, "Biomasa": b_flux, "Genes": ", ".join(cortes)})
    # Importante: Agregamos los cortes para prohibirlos en la siguiente copia
    acumulados.extend(cortes)

# Guardar y mostrar
df = pd.DataFrame(tier_list)
df.to_csv("tier_list_lactato_limpio.csv", index=False)
print("\n TIER LIST FINAL (MODO COPIA):")
print(df)