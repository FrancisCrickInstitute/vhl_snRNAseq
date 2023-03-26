# gene order file for inferCNV
devtools::load_all()

system(
  paste(
    "mkdir data/gencode ; cd data/gencode",
    "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz",
    "gunzip gencode.v43.basic.annotation.gtf.gz",
    "cat gencode.v43.basic.annotation.gtf | awk 'OFS=\"\t\" {if ($3==\"gene\") {print $1,$4-1,$5,$14}}' | tr -d '\";' | sort -k1,1 -k2,2n > gencode.v43.basic.annotation.bed",
    sep = ";"
  )
)
readr::read_tsv("data/gencode/gencode.v43.basic.annotation.bed",
                col_names = c("chr", "start", "end", "gene"),
                show_col_types = F) %>%
  # remove weird chromosomes / sex chromosomes
  dplyr::filter(!greplany(c("fix", "alt", "random", "chrUn", "chrM", "chrX", "chrY"), chr)) %>%
  # collapse to full range of each gene
  dplyr::group_by(gene) %>%
  dplyr::summarize(chr = unique(chr),
                   chr_n = as.numeric(gsub("chr", "", chr)),
                   start = min(start),
                   end = max(end)) %>%
  # remove duplicate genes
  dplyr::filter(dplyr::n() == 1) %>%
  # arrange by chr
  dplyr::arrange(chr_n) %>%
  dplyr::select(-chr_n) %>%
  readr::write_tsv("data/gencode/gencode.v43.basic.annotation_clean.bed", col_names = F)
