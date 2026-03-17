using TWWL
using Serialization
using DelimitedFiles

graphs = deserialize("mutag_graphs.jls")
H = 5
kernel = TwwlPrim(H)   # Usamos la versión que guarda las secuencias
fit!(kernel, graphs)
# Guardar las secuencias del primer grafo
seq_g1 = kernel.labelSeq[1]  # matriz (nodos × H+1)
writedlm("seq_julia_g1.csv", seq_g1, ',')
println("Secuencia guardada en seq_julia_g1.csv")