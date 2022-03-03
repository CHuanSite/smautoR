#' Structrally Masked Autoencoder
#'
#' R function of smautoPy for structrally masked autoencoder
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#'
#' @importFrom SJD
#'
#' @return A list containing the socre of each data set
#' @keywords
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#'               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#'               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#'               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
#' comp_num = c(2, 2, 2, 2, 2, 2, 2, 2, 2)
#'
#' res_sma = SMA(dataset, group, comp_num)
#'
#' @export

SMA <- function(dataset, group, comp_num, train_epoch = 1000){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    dataset = frameToMatrix(dataset)
    dataset = normalizeData(dataset)

    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])
    N_dataset = unlist(lapply(dataset, ncol))

    ## Output the component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : K){
        list_component[[i]] = matrix(0, nrow = p, ncol = comp_num[i])
        for(j in 1 : N){
            list_score[[j]][[i]] = matrix(0, nrow = comp_num[i], ncol = N_dataset[j])
        }
    }

    ## SMA called from Python by reticulate
    output_embedding = sma$smautoPy$StructuredMaskedAutoencoder(dataset, group, comp_num, N, p, c(p, 128, 32, 18), c(18, 32, 128, p), 100)

    ## compute the score for each dataset from SMA
    index = 1
    for(i in 1 : K){
        for(j in 1 : N){
            if (j %in% group[[i]]){
                list_score[[j]][[i]] = t(output_embedding[[j]][, index: (index + comp_num[i] - 1)])
            }
        }
        index = index + comp_num[i]
    }

    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = pveMultiple(dataset, group, comp_num, list_score, list_component)

    return(list(score_list = list_score))

}
