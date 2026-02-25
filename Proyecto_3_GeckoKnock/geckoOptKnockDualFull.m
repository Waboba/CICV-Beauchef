function [res] = geckoOptKnockDualFull(ecModel, biomassIdx, targetIdx, nKnockouts)
    %obtenemos dimensiones
    [nMets, nRxns] = size(ecModel.S);
    
    % 1. IDENTIFICAR CANDIDATAS
    todasLasProteinas = find(contains(ecModel.rxns, 'usage_prot_'));%buscamos todas las reacciones de uso de enzima
    excluir = contains(ecModel.rxns(todasLasProteinas), 'standard') | ...
              contains(ecModel.rxns(todasLasProteinas), 'draw_');%excluimos terminos que no son enzimas metabólicas reales
    candidatas = todasLasProteinas(~excluir);
    nCan = length(candidatas);
    
    % 2. BIG-M DINÁMICO
    M = max(ecModel.ub); 

    % 3. VARIABLES DEL PROBLEMA
    % [v (primal) | w (dual) | dLB | dUB | y (binarias)]
    nTotalVars = nRxns + nMets + nRxns + nRxns + nCan;
    v_idx   = 1:nRxns; %flujo de reacciones
    y_idx   = (nTotalVars - nCan + 1):nTotalVars;%variables binarias

    % 4. FUNCIÓN OBJETIVO (Bi-objetivo neutral)
    % Usamos un peso de 1 porque el "piso" de biomasa ya estará en el ecModel
    model.obj = zeros(nTotalVars, 1);
    model.obj(targetIdx) = 1.0;
    model.obj(biomassIdx) = 1.0; 
    model.modelsense = 'max';

    % 5. RESTRICCIONES
    % A. Primal: S*v = 0
    A1 = [ecModel.S, sparse(nMets, nTotalVars - nRxns)];
    rhs1 = zeros(nMets, 1);
    sense1 = repmat('=', nMets, 1);

    % B. Dual: S'*w + dLB - dUB = c (donde c_biom = 1)
    % Esto asegura que la célula siempre busque su óptimo interno
    c_vector = zeros(nRxns, 1);
    c_vector(biomassIdx) = 1;
    A2 = [sparse(nRxns, nRxns), ecModel.S', speye(nRxns), -speye(nRxns), sparse(nRxns, nCan)];
    rhs2 = c_vector;
    sense2 = repmat('=', nRxns, 1);

    % C. Acoplamiento de Knockouts: v_i <= M * y_i
    A3 = sparse(nCan, nTotalVars);
    for i = 1:nCan
        A3(i, candidatas(i)) = 1;
        A3(i, y_idx(i)) = -M;
    end
    rhs3 = zeros(nCan, 1);
    sense3 = repmat('<', nCan, 1);

    % D. Presupuesto de Knockouts
    A4 = [sparse(1, nTotalVars - nCan), ones(1, nCan)];
    rhs4 = nCan - nKnockouts;
    sense4 = '>';

    % Juntar Matriz
    model.A = [A1; A2; A3; A4];
    model.rhs = [rhs1; rhs2; rhs3; rhs4];
    model.sense = [sense1; sense2; sense3; sense4];

    % 6. BOUNDS Y TIPOS
    model.lb = -inf(nTotalVars, 1);
    model.lb(v_idx) = ecModel.lb; % Aquí vendrá el LB de biomasa del script principal
    model.lb((nRxns + nMets + 1):(nTotalVars - nCan)) = 0; 
    model.lb(y_idx) = 0;
    
    model.ub = inf(nTotalVars, 1);
    model.ub(v_idx) = ecModel.ub;
    model.ub(y_idx) = 1;

    model.vtype = [repmat('C', nTotalVars - nCan, 1); repmat('B', nCan, 1)];

    % 7. CONFIGURACIÓN DEL SOLVER
    params.Presolve = 0;   
    params.TimeLimit = 600;
    params.MIPFocus = 3;   

    res = gurobi(model, params);
end