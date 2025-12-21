
library(ELVAR)
library(Seurat)
library(Matrix)
library(dplyr)
library(igraph)
library(MASS)
library(ggplot2)


run_da_elvar <- function(obj,
                      resolution,
                      threshold,
                      alpha,
                      Vattr.name) {
    
    # ind the names of all graphs ending with "_nn"
    nn_graphs <- grep("_nn$", names(obj@graphs), value = TRUE)
     if (length(nn_graphs) == 0) {
        stop("No *_nn graph found in obj@graphs")
    }
    
    # Select the adjacency matrix priority
    if ("integrated_nn" %in% nn_graphs) {
        chosen_graph <- "integrated_nn"
    } else if ("RNA_nn" %in% nn_graphs) {
        chosen_graph <- "RNA_nn"
    } else if (length(nn_graphs) == 1) {
        chosen_graph <- nn_graphs[1]
    } else {
        # If there are multiple but no integrated_nn and RNA_nn, select the first one
        chosen_graph <- nn_graphs[1]
        warning(paste("Multiple *_nn graphs found, using", chosen_graph))
    }
    
    message("Using graph: ", chosen_graph)

    adj.m <- as.matrix(obj@graphs[[chosen_graph]])

    diag(adj.m) <- 0

    gr.o <- graph.adjacency(adj.m, mode="undirected")

    vertex_attr(gr.o, name="condition") <- obj@meta.data$condition

    print(is.connected(gr.o))

    eva.o <- Eva_partitions(gr.o, resolution=resolution, threshold=threshold,
                            alpha=alpha, Vattr.name=Vattr.name)

    print(table(eva.o$CommunityMembers, obj@meta.data$condition))

    sigcl.o <- ProcessEVA(eva.o, obj, attrName="condition")

    print(sigcl.o$pvalEnr)
    print(sigcl.o$sigClust)

    print(table(obj@meta.data$condition[unlist(sigcl.o$cellsMrg)],
                obj@meta.data$orig.ident[unlist(sigcl.o$cellsMrg)]))


    da.o <- DoDA(obj, sigcl.o, varDA="condition", varREP="sample")
    print(da.o$stat)


  # Return the relevant result you want to store
    return(list(sigcl.o = sigcl.o, da_o = da.o))


}




