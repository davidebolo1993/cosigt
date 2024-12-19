library(data.table)
library(dbscan)
library(rjson)
library(reshape2)
library(reshape2)
library(NbClust)

input_file<-"filt.tsv.gz"
df<-fread(input_file)

for (d in c("euclidean.dist","jaccard.dist","cosine.dissim","manhattan.dist")) {

    regularMatrix <- acast(df, group.a ~ group.b, value.var = d)
    distanceMatrix<-as.dist(regularMatrix)
    pdf(paste0("knn.",d,".pdf"))
    kNNdistplot(distanceMatrix,k=2)
    dev.off()
    kNN_distances <- kNNdist(distanceMatrix, k = 2)
    sorted_kNN <- sort(kNN_distances)
    first_derivative <- diff(sorted_kNN)
    # Step 2: Compute the second derivative
    second_derivative <- diff(first_derivative)
    # Step 3: Identify the index with the maximum second derivative
    optimal_index <- which.max(second_derivative)
    # Step 4: Retrieve the corresponding `eps` value
    optimal_eps <- sorted_kNN[optimal_index + 1]  # +1 d
    db<-dbscan(distanceMatrix,minPts=3, eps=4.3)
    cl<-db$cluster
    names(cl)<-labels(distanceMatrix)
    res.list <- lapply(split(cl, names(cl)), unname)
    named_res <- lapply(cl, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
    jout <- toJSON(named_res)
    # Write JSON output
    output_file<-paste0("dbscan.",d,".json")
    write(jout, output_file)


    max_cluster <- round(length(unique(df$group.a)) / 5) ##control
    res <- NbClust(diss = distanceMatrix, method = "average", index = "silhouette", 
                distance = NULL, max.nc = max_cluster)$Best.partition

    # Format results
    res.list <- lapply(split(res, names(res)), unname)
    named_res <- lapply(res.list, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
    jout <- toJSON(named_res)

    # Write JSON output
    output_file<-paste0("agglomerative.",d,".json")
    write(jout, output_file)

}