library(fagin2)

get_arabidopsis_config <- function(){

  # create a default configuration
  con <- config()

  # set parameter
  con@synder@offsets = c(1L,1L) # offsets for mummer, c(0L,1L) for satsuma
  con@synder@trans = "p" # percent identity transform (0-100), or "d" for proportion transform (0-1)
  # alnrate is the proportion of target sequence match to query sequence for alignment
  # default cutoff for alnrate are 0
  con@alignment@alnrate@prot2prot=0.5 # protein vs protein
  con@alignment@alnrate@prot2allorf=0.5 # protein vs all ORF
  con@alignment@alnrate@prot2transorf=0.5 #protein vs ORF on the mRNA
  con@alignment@alnrate@dna2dna=0.5 # DNA vs DNA

  # focal species can be one or more species (max: the total number of species in tree file)
  con@input@focal_species = c("Arabidopsis_thaliana")

  # set path for input files
  # name of list should be same as the species in tree
  con@input@gff <- list(
    "Arabidopsis_thaliana"   = "Arabidopsis_thaliana_annotation.gff3"
    , "Arabidopsis_lyrata"   = "Arabidopsis_lyrata_annotation.gff3"
    , "Arabidopsis_halleri" = "Arabidopsis_halleri_annotation.gff3"
    , "Eutrema_salsugineum" = "Eutrema_salsugineum_annotation.gff3"
    , "Brassica_napus" = "Brassica_napus_annotation.gff3"
    , "Brassica_olerecea" = "Brassica_olerecea_annotation.gff3"
    , "Athal_An-1" = "Athal_An-1_annotation.gff3"
    , "Athal_C24" = "Athal_C24_annotation.gff3"
    , "Athal_Cvi" = "Athal_Cvi_annotation.gff3"
    , "Athal_Eri1" = "Athal_Eri1_annotation.gff3"
    , "Athal_Kyo" = "Athal_Kyo_annotation.gff3"
    , "Athal_Ler" = "Athal_Ler_annotation.gff3"
    , "Athal_Sha" = "Athal_Sha_annotation.gff3"
  )
  con@input@fna <- list(
    "Arabidopsis_thaliana"   = "Arabidopsis_thaliana_genome.fna"
    , "Arabidopsis_lyrata"   = "Arabidopsis_lyrata_genome.fna"
    , "Arabidopsis_halleri" = "Arabidopsis_halleri_genome.fna"
    , "Eutrema_salsugineum" = "Eutrema_salsugineum_genome.fna"
    , "Brassica_napus" = "Brassica_napus_genome.fna"
    , "Brassica_olerecea" = "Brassica_olerecea_genome.fna"
    , "Athal_An-1" = "Athal_An-1_genome.fna"
    , "Athal_C24" = "Athal_C24_genome.fna"
    , "Athal_Cvi" = "Athal_Cvi_genome.fna"
    , "Athal_Eri1" = "Athal_Eri1_genome.fna"
    , "Athal_Kyo" = "Athal_Kyo_genome.fna"
    , "Athal_Ler" = "Athal_Ler_genome.fna"
    , "Athal_Sha" = "Athal_Sha_genome.fna"
  )

  # syn path is a two-dimensional list
  # first level: focal species; second level: target species
  con@input@syn <- list(
    Arabidopsis_thaliana = list(
      "Arabidopsis_lyrata"   = "Arabidopsis_thaliana_Arabidopsis_lyrata.syn"
      , "Arabidopsis_halleri" = "Arabidopsis_thaliana_Arabidopsis_halleri.syn"
      , "Eutrema_salsugineum" = "Arabidopsis_thaliana_Eutrema_salsugineum.syn"
      , "Brassica_napus" = "Arabidopsis_thaliana_Brassica_napus.syn"
      , "Brassica_olerecea" = "Arabidopsis_thaliana_Brassica_olerecea.syn"
      , "Athal_An-1" = "Arabidopsis_thaliana_Athal_An-1.syn"
      , "Athal_C24" = "Arabidopsis_thaliana_Athal_C24.syn"
      , "Athal_Cvi" = "Arabidopsis_thaliana_Athal_Cvi.syn"
      , "Athal_Eri1" = "Arabidopsis_thaliana_Athal_Eri1.syn"
      , "Athal_Kyo" = "Arabidopsis_thaliana_Athal_Kyo.syn"
      , "Athal_Ler" = "Arabidopsis_thaliana_Athal_Ler.syn"
      , "Athal_Sha" = "Arabidopsis_thaliana_Athal_Sha.syn"
    )
  )
  con@input@tree <- "tree"
  con@input@gene_list <- list(
    "Arabidopsis_thaliana" = "Arabidopsis_thaliana_annotation_gene_list.txt"
    , "Brassica_napus"   =  "Brassica_napus_annotation_gene_list.txt"
    , "Arabidopsis_halleri" = "Arabidopsis_halleri_annotation_gene_list.txt"
    , "Eutrema_salsugineum" = "Eutrema_salsugineum_annotation_gene_list.txt"
    , "Arabidopsis_lyrata" = "Arabidopsis_lyrata_annotation_gene_list.txt"
    , "Brassica_olerecea" = "Brassica_olerecea_annotation_gene_list.txt"
    , "Athal_An-1" = "Athal_An-1_annotation_gene_list.txt"
    , "Athal_C24" = "Athal_C24_annotation_gene_list.txt"
    , "Athal_Cvi" = "Athal_Cvi_annotation_gene_list.txt"
    , "Athal_Eri1" = "Athal_Eri1_annotation_gene_list.txt"
    , "Athal_Kyo" = "Athal_Kyo_annotation_gene_list.txt"
    , "Athal_Ler" = "Athal_Ler_annotation_gene_list.txt"
    , "Athal_Sha" = "Athal_Sha_annotation_gene_list.txt"
  )
  con@archive = "arabidopsis_thaliana_fagin_archive" # dir for output files

  validate_config(con)
  con
}

con <- get_arabidopsis_config() # create revised configuration
m <- run_fagin_parallel(con, cores=16, cl.type="FORK")
