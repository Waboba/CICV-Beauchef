using TWWL
using Graphs
using MetaGraphsNext
using DelimitedFiles
using Serialization

# Descargar dataset MUTAG
inp = downloadDataset("MUTAG")
graphs = inp[:data]  # lista de MetaGraph

# Crear carpeta de salida (opcional, para los archivos de texto)
output_dir = "mutag_graphs"
mkpath(output_dir)

# Para cada grafo, exportar a un archivo de texto (por si quieres inspeccionarlos)
for (i, mg) in enumerate(graphs)
    n = nv(mg)
    labels = [mg[label_for(mg, v)][:label] for v in 1:n]
    edge_list = [(e.src, e.dst) for e in edges(mg.graph)]
    
    filename = joinpath(output_dir, "graph_$(i).txt")
    open(filename, "w") do f
        println(f, n)
        println(f, join(labels, " "))
        for (u, v) in edge_list
            println(f, "$u $v")
        end
    end
end

# Guardar los grafos originales en un archivo serializado
serialize("mutag_graphs.jls", graphs)
println("Exportación completada. Grafos originales guardados en mutag_graphs.jls")