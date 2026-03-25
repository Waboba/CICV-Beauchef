

# %%
def ENM(modelo_base, mu_medido, lista_proteinas, L=None, medio=None):
    # Determinamos qué nutrientes evaluar: las llaves del medio o la lista L
    nutrientes_a_evaluar = list(medio.keys()) if medio is not None else L
    if nutrientes_a_evaluar is None:
        raise ValueError("Debes proporcionar al menos una lista L o un diccionario de medio.")
    
    bool_array = np.zeros(len(nutrientes_a_evaluar), dtype=int)
    minimum_uptake = np.zeros(len(nutrientes_a_evaluar))

    # Creamos el modelo base para esta corrida y aplicamos el medio
    model1 = modelo_base.copy()
    if medio is not None:
        model1.medium = medio

    sol_opt = model1.optimize()
    

 
    # Buscamos la ID de la biomasa
    biomasa = modelo_base.reactions.get_by_id(modelo_base.objective.expression.args[0].args[1].name)
    id_biomasa= biomasa.id
    p_medida = list(lista_proteinas.values())

    # Usamos la lista de nutrientes detectada al inicio
    for i, nutriente_objetivo in enumerate(nutrientes_a_evaluar):
        with model1 as modelo_temp:



            for nutriente in nutrientes_a_evaluar:
                rxn_1 = modelo_temp.reactions.get_by_id(nutriente)
                rxn_1.bounds=-1000,1000
            
            # v_biom = mu_obs
            rxn_biom = modelo_temp.reactions.get_by_id(id_biomasa)
            rxn_biom.bounds = (mu_medido, mu_medido)

            # v_prot = p_obs 
            for j, rid in enumerate(lista_proteinas):
                rxn = modelo_temp.reactions.get_by_id(rid)
                rxn.bounds = (p_medida[j], p_medida[j])


            


           # Pqueño arreglo para liealizar el valor absoluto
            rxn_2 = modelo_temp.reactions.get_by_id(nutriente_objetivo) 
            a = modelo_temp.solver.interface.Variable(f"abs_var_{nutriente_objetivo}", lb=0)
            modelo_temp.solver.add(a)

            modelo_temp.solver.add(modelo_temp.solver.interface.Constraint(
             a - rxn_2.flux_expression, lb=0, ub=None, name=f"abs_ge_{nutriente_objetivo}"  ))

            modelo_temp.solver.add(modelo_temp.solver.interface.Constraint(
            a + rxn_2.flux_expression, lb=0, ub=None, name=f"abs_le_{nutriente_objetivo}" ))

            modelo_temp.objective = a
            modelo_temp.objective.direction = 'min'    

            solucion = modelo_temp.optimize()


            
            if solucion.status == 'optimal':
                # Usamos una tolerancia pequeña para el cero numérico
                if abs(solucion.objective_value) < 1e-8:
                    bool_array[i] = 1 # No esencial
                
                minimum_uptake[i] = solucion.fluxes[nutriente_objetivo]
            else:
                minimum_uptake[i] = np.nan

    return bool_array, minimum_uptake




# %%
def UOF(modelo_base, mu_medido, lista_proteinas, matriz, valores_minimos, L=None, medio=None):


    #Obtenemos los utrientes del diccionario o de la lista
    nutrientes_a_evaluar = list(medio.keys()) if medio is not None else L

   
    #Copiamos el modelo el insertamos el meio en caso de haber recibido uno
    model2 = modelo_base.copy()
    if medio is not None:
        model2.medium = medio

    sol_opt = model2.optimize()
    biom_id = model2.reactions.get_by_id(model2.objective.expression.args[0].args[1].name).id 
    p_medida = list(lista_proteinas.values())


   
    nutrientes = np.array(nutrientes_a_evaluar)
    esenciales = nutrientes[matriz == 0]
    no_esenciales = nutrientes[matriz == 1]
    minimos_esenciales = valores_minimos[matriz == 0]

    soluciones = []

    for i, nut_no_esencial in enumerate(no_esenciales):
        # 4. Uso de 'with' para velocidad
        with model2 as modelo_temp2:
            
            # Fijar Biomasa usando el ID de texto
            rxn_biom = modelo_temp2.reactions.get_by_id(biom_id)
            rxn_biom.bounds = (mu_medido, mu_medido)

            # Fijar Proteínas/Productos
            for j, rid in enumerate(lista_proteinas):
                rxn = modelo_temp2.reactions.get_by_id(rid)
                valor = p_medida[j]
                rxn.bounds = (valor, valor)

            # Restricciones de esenciales
            
            for k, nombre_es in enumerate(esenciales):
                rx = modelo_temp2.reactions.get_by_id(nombre_es)
                if abs(rx.upper_bound)<abs(minimos_esenciales[k]):
                    print("it happened")
                    rx.bounds = minimos_esenciales[k],minimos_esenciales[k]

            # Configuración del nutriente no esencial objetivo
            rxn_act = modelo_temp2.reactions.get_by_id(nut_no_esencial)
            rxn_act.bounds = (-1000, 1000)

            # Valor absoluto linealizado
            a = modelo_temp2.solver.interface.Variable(f"abs_uof_{nut_no_esencial}", lb=0)
            modelo_temp2.solver.add(a)
            
            modelo_temp2.solver.add(modelo_temp2.solver.interface.Constraint(
                a - rxn_act.flux_expression, lb=0, name=f"uof_ge_{nut_no_esencial}"))
            modelo_temp2.solver.add(modelo_temp2.solver.interface.Constraint(
                a + rxn_act.flux_expression, lb=0, name=f"uof_le_{nut_no_esencial}"))

            modelo_temp2.objective = a
            modelo_temp2.objective.direction = 'min'
            
            # Optimización
            opt = modelo_temp2.optimize()
            soluciones.append((nut_no_esencial, opt))

    return esenciales, soluciones


# %%
def mapa_calor(modelo_base, mu_medido, lista_proteinas, L=None, medio=None):
    # Llamada a ENM para clasificar nutrientes
    bool_array, minimun_uptakes = ENM(modelo_base, mu_medido, lista_proteinas, L=L, medio=medio)
    print(bool_array)
    print(minimun_uptakes)
    
    # Llamada a UOF usando los resultados de ENM
    esenciales_nombres, soluciones = UOF(modelo_base, mu_medido, lista_proteinas, bool_array, minimun_uptakes, L=L, medio=medio)

    print(esenciales_nombres)
    print(soluciones)

    # Separación de nutrientes para la matriz (basado en bool_array)
    nutrientes_arr = np.array(list(medio.keys()) if medio is not None else L)
    esenciales = nutrientes_arr 
    no_esenciales = nutrientes_arr[bool_array == 1]

    # Crear DataFrame para Precios Sombra
    matriz_dual = pd.DataFrame(index=no_esenciales, columns=esenciales)


    for nut_no_es, sol in soluciones:
        if sol.status == 'optimal':
            precios_sombra = sol.shadow_prices
            for nut_es in esenciales:
                reaccion = modelo_base.reactions.get_by_id(nut_es)
                met_id = reaccion.reactants[0].id
                matriz_dual.loc[nut_no_es, nut_es] = abs(precios_sombra[met_id])
        else:
            matriz_dual.loc[nut_no_es, :] = np.nan

    matriz_plot = matriz_dual.astype(float)

    plt.figure(figsize=(12, 8))
    sns.heatmap(matriz_plot, annot=True, cmap='YlOrRd', fmt=".2f", linewidths=.5)
    plt.title('Mapa de Calor: Precios Sombra')
    plt.xlabel('Nutrientes Esenciales')
    plt.ylabel('Nutrientes No Esenciales')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()



