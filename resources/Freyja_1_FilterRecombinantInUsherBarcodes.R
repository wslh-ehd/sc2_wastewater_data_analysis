rm(list=ls(all=TRUE))
ls()
set.seed(123)



require(tidyverse)
require(strex)
require(data.table)
require(stringi)
require(jsonlite)
require(yaml)
'%!in%' <- function(x,y)!('%in%'(x,y)) 
'%!like%' <- function(x,y)!('%like%'(x,y)) 






#################################################################################
#### Import data
#################################################################################

# Import and reformat alias_key.json
alias_key <- jsonlite::read_json("alias_key.json", simplifyVector = TRUE)
alias_key <- cbind(names(alias_key), as.data.frame(matrix(c(alias_key, rep(NA, length(alias_key) %% length(alias_key))), length(alias_key)))); names(alias_key)<-c("original", "temp")

alias_key<-as.data.frame(lapply(alias_key, gsub, pattern='"', replacement=""))
alias_key<-as.data.frame(lapply(alias_key, gsub, pattern="[*]", replacement=""))
alias_key<-as.data.frame(lapply(alias_key, gsub, pattern=" ", replacement=""))

alias_key <- alias_key %>%
  as.data.frame() %>%
  dplyr::mutate(n = stringr::str_count(as.character(temp), ",") + 1 ,
                final = gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", temp, perl=T),
                final = ifelse(final == "", temp, final),
                final_simplified = ifelse(n==1, sub("\\..*", "", final), "")) %>%
  dplyr::filter(original %!in% c("A", "B")) %>% # Remove A and B from alias_keys - Important for Freyja_3_Prep_visual.R
  dplyr::select(original, final, final_simplified, n)



# Import and reformat clade_display_names.yml 
clade_display_names <- yaml::read_yaml("clade_display_names.yml")
clade_display_names<-do.call("rbind", clade_display_names) %>%
  data.frame() %>% # Convert into df
  tidyr::separate(".", into=c("Nextstrain", "Pango"), sep=" ", convert = TRUE, fill = "right") %>%
  dplyr::mutate(Pango = gsub('[)]', "", gsub('[(]', "", Pango)),
                Pango_simplified = sub("\\..*", "", Pango))



# Import usher_barcodes
usher_barcodes<-read.delim("usher_barcodes.csv", header = FALSE)







#################################################################################
#### Identify Nextstrain recombinant
#################################################################################
nextstrain.recombinant<-clade_display_names %>%
  dplyr::filter(stringr::str_detect(Pango_simplified, "^X")) %>%
  dplyr::select(Pango_simplified) %>%
  unique()



#################################################################################
#### Identify recombinants
#################################################################################

#### Identify recombinants ####
alias_key_list.recombinant.1 <- alias_key %>%
  dplyr::mutate(n = as.numeric(as.character(n))) %>%
  dplyr::filter(n > 1) %>%
  unique()

#### Identify recombinants that do not start with X ####
alias_key_list.recombinant.2 <- alias_key %>%
  dplyr::mutate(n = as.numeric(as.character(n))) %>%
  dplyr::filter(n == 1, 
                stringr::str_detect(final, "^X")) %>%
  unique()


#### Merge 2x datasets ####
alias_key_list.recombinant <- rbind(alias_key_list.recombinant.1, alias_key_list.recombinant.2)



#################################################################################
#### Identify recombinants to discard do not start with "^X"
#################################################################################

list.recombinant.to.discard<-alias_key_list.recombinant %>%
  dplyr::filter(original %!in% nextstrain.recombinant$Pango_simplified) %>% # Remove recombinants (most likely starting by "X")
  dplyr::filter(final_simplified %!in% nextstrain.recombinant$Pango_simplified) %>% # Remove the recombinants that do not start with "X"
  dplyr::select(original) %>%
  unique()



#################################################################################
#### USHER DATABASE: Filter out the recombinants (not described as Nextstrain clades) from usher_barcodes
#################################################################################

# Reformat usher_barcodes
usher_barcodes.processed<-stringr::str_split_fixed(usher_barcodes$V1, pattern = ',', n = 2) %>%
  as.data.frame() %>%
  dplyr::mutate(simplified = sub("\\..*", "", V1))

# Remove non Nextstrain recombinant clades
usher_barcodes.filtered<-usher_barcodes.processed %>%
  dplyr::filter(simplified %!in% list.recombinant.to.discard$original) %>%
  dplyr::select(-simplified)






# Export filtered usher_barcodes
write.table(usher_barcodes.filtered, "usher_barcodes_FilteredRecombinants.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(usher_barcodes.filtered$V1, "usher_barcodes_FilteredRecombinants_listVariants.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(usher_barcodes.processed$V1, "usher_barcodes_listVariants.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(clade_display_names, "clade_display_names.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(alias_key, "alias_key.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(alias_key_list.recombinant, "alias_key_recombinant.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
