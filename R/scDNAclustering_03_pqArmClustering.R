##### 03: pqArm Clustering

# 3.1: CN_template() build the CN_seq template, stands the range for each row in CN_seq
CN_template <- function(input, pqArm_file){
  CN_tem <- data.frame(seqnames(input[[1]]$bins), ranges(input[[1]]$bins))
  CN_tem <- CN_tem %>%
    setNames(c("chr", "start", "end", "width"))

  # Add p/q arm information on the template
  # 非p及q
  pqArm_range <- pqArm_file.remake(FILE = pqArm_file)
  pqArm_range <- pqArm_range %>%
    pqArm_file.pq(Template = .) %>%
    filter(arm == "p")

  CN_tem_pq <- NULL
  for (i in 1:nrow(pqArm_range)){
    CN_tem_pq <- CN_tem %>%
      filter(chr %in% pqArm_range$chr[i],
             start >= pqArm_range$start[i],
             end <= pqArm_range$end[i]) %>%
      mutate(arm = pqArm_range$arm[i]) %>%
      rbind(CN_tem_pq)
  }

  CN_tem_pq <- left_join(CN_tem, CN_tem_pq, by = c("chr", "start", "end", "width"))
  CN_tem_pq <- CN_tem_pq %>%
    replace_na(list(arm = "q")) %>%
    arrange(chr, start)

  return(CN_tem_pq)
}


# 3.1.1: pqArm_file.remake() remake file format from UCSC
# FILE = "~/its00/DNA_seq/scDNA_clustering/Template_data/hg38_cytoBand.txt.gz"
pqArm_file.remake <- function(FILE){
  x <- read.table(gzfile(FILE),sep="\t", col.names = c("chr", "start", "end", "name","gieStain"))
  x <- x %>%
    filter(!grepl("_", chr))

  # set chr levels
  vec <- unique(x$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  Levels <- c(nums, vec[!grepl("\\d", vec)])


  x <- x %>%
    mutate(chr = factor(chr, levels = Levels),
           start = start + 1,
           arm  = substring(name, 1, 1),
           arm_category = paste0(chr, arm))
  Sum_x <- x %>%
    group_by(arm_category) %>%
    slice_head(n = 1) %>%
    as.data.frame()
  Sum_x <- x %>%
    group_by(arm_category) %>%
    slice_tail(n = 1) %>%
    as.data.frame() %>%
    rbind(Sum_x) %>%
    arrange(chr = factor(chr, levels = Levels), start)

  return(Sum_x)
}


# 3.1.2: pqArm_file.pq() remake pqarm template into new format
pqArm_file.pq <- function(Template){
  x <- Template %>%
    setNames(c("chr", "ChromStart", "ChromEnd", "name", "gieStain", "arm", "arm_category")) %>%
    group_by(arm_category) %>%
    summarise(start = min(ChromStart),
              end = max(ChromEnd)) %>%
    as.data.frame()
  x <- x %>%
    mutate(arm = str_sub(arm_category, -1),
           chr = str_sub(arm_category, end = -2)) %>%
    select(c("chr", "start", "end", "arm"))

  return(x)
}


# 3.1.3: pqArm_file.cen() remake centromere template into new format (X)
pqArm_file.cen <- function(FILE){
  x <- read.table(gzfile(FILE),sep="\t", col.names = c("chr", "ChromStart", "ChromEnd", "name","gieStain"))
  x <- x %>% filter(!grepl("_", chr))

  # set chr levels
  vec <- unique(x$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  Levels <- c(nums, vec[!grepl("\\d", vec)])

  x <- x %>%
    mutate(cen_category = paste0(chr, gieStain)) %>%
    group_by(cen_category) %>%
    summarise(Start = min(ChromStart),
              End  = max(ChromEnd)) %>%
    as.data.frame()

  x <- x %>%
    filter(cen_category %in% paste0(rep(Levels,each = 2), c("acen","gvar"))) %>%
    mutate(cen = str_sub(cen_category, -4),
           chr = str_sub(cen_category, end = -5)) %>%
    group_by(chr) %>%
    summarise(MaskStart = min(Start),
              MaskEnd  = max(End)) %>%
    arrange(chr = factor(chr, levels = Levels), MaskStart) %>%
    select(c("chr", "MaskStart", "MaskEnd"))

  return(x)
}


# 3.2: pqArm_CN() transform CN matrix to Del/Neu/Amp and based on Arm level to smooth the CN
pqArm_CN <- function(input, Cluster_label, Clustering_output, pqArm_file){
  selected_files <- subset(Clustering_output, Clustering_output$cluster%in%c(Cluster_label))

  CN_matrix_temp <- CN_template(input = input, pqArm_file = pqArm_file)
  CN_matrix <- CN_seq(input = input, Template = selected_files$cell)
  CN_matrix <- pqArm_DelNeuAmp(matrix = CN_matrix)  # Transform CN to 1, 2, 3

  CN_binsLevel <- CN_matrix_temp %>%
    group_by(chr, arm) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    as.data.frame()
  CN_binsLevel$start <- sapply(1:nrow(CN_binsLevel), function(x){
    sum(CN_binsLevel$Freq[1:x-1])+1
  })
  CN_binsLevel$end <- sapply(1:nrow(CN_binsLevel), function(x){
    sum(CN_binsLevel$Freq[1:x])
  })

  Smooth_pqCN <- NULL
  for(i in 1:nrow(CN_binsLevel)){
    pq_CNmatrix <- CN_matrix[CN_binsLevel$start[i]:CN_binsLevel$end[i], ]
    pq_CN <- sapply(1:ncol(pq_CNmatrix), function(x){
      freq <- table(pq_CNmatrix[ , x]) %>%
        as.data.frame() %>%
        arrange(desc(Freq)) %>%
        pull(Var1) %>%
        as.character() %>%
        as.integer()
      first_element <- freq[1]
    })

    Smooth_pqCN <- Smooth_pqCN %>%
      rbind(pq_CN)
  }

  Smooth_pqCN <- Smooth_pqCN %>%
    as.data.frame() %>%
    setNames(c(colnames(CN_matrix))) %>%
    `rownames<-`(paste0(rep(levels(CN_matrix_temp$chr), each = 2), rep(c("p", "q"), 11)))

  return(Smooth_pqCN)
}


# 3.2.1: pqArm_DelNeuAmp() transform CN matrix to Del/Neu/Amp three types
pqArm_DelNeuAmp <- function(matrix){
  new_matrix <- matrix(0, nrow(matrix), ncol(matrix))

  new_matrix[matrix < 2] <- 1
  new_matrix[matrix == 2] <- 2
  new_matrix[matrix > 2] <- 3

  colnames(new_matrix) <- colnames(matrix)


  return(new_matrix)
}


# 3.3: pqArm_clustering() merge each row as vector to seperate the clsuters
pqArm_clustering <- function(matrix, Label){
  cluster <- sapply(1:ncol(matrix), function(x){
    paste(matrix[ , x], collapse = "_")
  })

  cluster <- cluster %>%
    data.frame(pqArm_pattern = .) %>%
    cbind(cellID = colnames(matrix)) %>%
    mutate(cluster = Label)


  return(cluster)
}


# 3.4: pqArm_clustering_summary() return pqArm clustering output
pqArm_clustering_summary <- function(matrix, Label){
  cluster_table <- table(matrix$pqArm_pattern) %>%
    as.data.frame() %>%
    arrange(desc(Freq)) %>%
    mutate(cluster = Label,
           pqArm_cluster = c(1:nrow(.))) %>%
    # filter(Freq > 1) %>%
    setNames(c("pqArm_pattern", "pqArm_pattern_cellnum", "cluster", "pqArm_cluster"))

  return(cluster_table)
}
