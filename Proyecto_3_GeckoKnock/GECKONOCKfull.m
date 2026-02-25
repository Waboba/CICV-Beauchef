
% Script para la identificación de knockouts enzimáticos mediante OptKnock Dual
% aplicado a modelos GECKO (Enzyme-Constrained).
% Objetivo: Maximizar la producción acoplada al crecimiento.

clear; clc;

% --- 1. CONFIGURACIÓN E INICIALIZACIÓN ---
fprintf('Cargando adaptador y modelo GECKO...\n');
adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params_g = ModelAdapter.getParameters();

% Carga del modelo (Asegurarse de que el archivo .yml esté en tu path)
ecModel = loadEcModel('ecYeastGEM.yml'); 

% --- 2. CONFIGURACIÓN DE ÍNDICES Y TARGET ---
biomassIdx = find(strcmp(ecModel.rxns, params_g.bioRxn)); 

% Define el ID de la reacción que se quiera producir:
% r_1546: L-Lactato | r_1914: L-Lisina | r_2056: Succinato
targetID = 'r_1914'; 
targetIdx = find(strcmp(ecModel.rxns, targetID));

if isempty(targetIdx)
    error('El ID del target no fue encontrado en el modelo.');
end

% --- 3. RESTRICCIONES AMBIENTALES (Condiciones de Cultivo) ---
% Glucosa (r_1714): -10 mmol/gDWh
ecModel.lb(find(strcmp(ecModel.rxns, 'r_1714'))) = -10; 
% Oxígeno (r_1992): -20 mmol/gDWh (Condiciones aeróbicas)
ecModel.lb(find(strcmp(ecModel.rxns, 'r_1992'))) = -20; 

% --- 4. FIJAR PISO DE BIOMASA (Garantía de Viabilidad) ---
% Calculamos el máximo teórico para definir un 50% de crecimiento mínimo.
% Esto asegura que el diseño sea realista.
fprintf('Calculando biomasa máxima teórica...\n');
solFBA = solveLP(ecModel);
biomasa_minima = solFBA.f * 0.50;%variable a ajustar dependiendo del crecimiento mínimo que se necesite.
ecModel.lb(biomassIdx) = biomasa_minima; 

fprintf('Configuración final:\n');
fprintf(' - Target: %s\n', targetID);
fprintf(' - Crecimiento mínimo forzado: %.4f h^-1 (50%% del max)\n', biomasa_minima);

% --- 5. EJECUCIÓN DE OPTKNOCK DUAL ---
nKnockouts = 3; % Número máximo de enzimas a eliminar
fprintf('\nIniciando OptKnock Dual Completo (Problema Bi-nivel)...\n');
fprintf('Este proceso puede tomar varios minutos debido a la complejidad de la matriz...\n');

res = geckoOptKnockDualFull(ecModel, biomassIdx, targetIdx, nKnockouts);

% --- 6. BLOQUE DE REPORTE DE RESULTADOS ---
if isstruct(res) && isfield(res, 'status') && (strcmp(res.status, 'OPTIMAL') || strcmp(res.status, 'SUBOPTIMAL'))
    
    % Recuperar flujos reales del vector de variables x
    flujo_biomasa = res.x(biomassIdx);
    flujo_target = res.x(targetIdx);
    
    fprintf('\n========================================\n');
    fprintf('   RESULTADOS DEL DISEÑO OPTKNOCK\n');
    fprintf('========================================\n');
    fprintf('Producción de %s: %.6f mmol/gDWh\n', targetID, flujo_target);
    fprintf('Crecimiento celular: %.6f h^-1\n', flujo_biomasa);
    fprintf('Estatus del Solver: %s\n', res.status);
    
    % Identificar qué enzimas fueron "apagadas" (binarias y_i < 0.5)
    todasLasProteinas = find(contains(ecModel.rxns, 'usage_prot_'));
    excluir = contains(ecModel.rxns(todasLasProteinas), 'standard') | ...
              contains(ecModel.rxns(todasLasProteinas), 'draw_');
    candidatas = todasLasProteinas(~excluir);
    nCan = length(candidatas);
    
    % Las variables binarias están al final del vector x
    y_values = res.x(end - nCan + 1 : end);
    cortadas_idx = candidatas(y_values < 0.5);
    
    if isempty(cortadas_idx)
        fprintf('\nResultado: No se identificaron knockouts beneficiosos.\n');
    else
        fprintf('\nDISEÑO DE ENZIMAS PARA ELIMINACIÓN:\n');
        for i = 1:length(cortadas_idx)
            idx = cortadas_idx(i);
            geneName = strrep(ecModel.rxnNames{idx}, 'usage_prot_', '');
            fprintf(' [%d] ID: %s | Enzima/Proteína: %s\n', i, ecModel.rxns{idx}, geneName);
        end
        fprintf('\nNota: Este diseño garantiza que la producción es obligatoria\n');
        fprintf('si la célula crece a la velocidad indicada.\n');
    end
    
else
    fprintf('\nEl solver no pudo encontrar una solución óptima.\n');
    if isfield(res, 'status')
        fprintf('Estatus: %s\n', res.status);
    end
end