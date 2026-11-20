# ==============================================================================
# 00_data_audit.R -- integrity checks on the raw data (run first)
# Outputs: results/tables/audit_*.csv
# ==============================================================================
source("config.R"); source("R/functions.R")
raw <- read.table(CFG$data_txt, header = TRUE, sep = "\t", check.names = FALSE)
raw$plot <- as.integer(sub("^P(\\d+)S\\d+$", "\\1", raw$Unidad))

# 1. leaves per campaign (raw file)
a1 <- as.data.frame(table(campaign = raw$Dia)); names(a1)[2] <- "n_raw"
# 2. duplicated keys
key <- paste(raw$Dia, raw$Unidad)
a2 <- raw[key %in% key[duplicated(key)], c("Dia", "Unidad", "Bloque_Num", "Tratamiento", "Chla", "Chlb")]
a2$row_in_file <- which(key %in% key[duplicated(key)])
# 3. leaves that are missing in a full 20 plots x 4 leaves design
full <- paste(rep(names(DDC_MAP), each = 80), paste0("P", rep(1:20, each = 4), "S", 1:4))
a3 <- data.frame(missing_key = setdiff(full, key))
# 4. block x treatment table (a balanced RCBD has exactly one plot per cell)
pl <- unique(raw[, c("plot", "Bloque_Num", "Tratamiento")])
a4 <- as.data.frame.matrix(table(block = pl$Bloque_Num, treatment = pl$Tratamiento))
a4 <- cbind(block = rownames(a4), a4)
# 5. the observation the original script dropped without explanation (row 468)
x4 <- raw[raw$Dia == "Dia 4", ]; z <- (x4$Chla - mean(x4$Chla)) / sd(x4$Chla)
a5 <- data.frame(key = "Dia 4 P6S1", Chla = x4$Chla[x4$Unidad == "P6S1"], z_score_within_campaign = round(z[x4$Unidad == "P6S1"], 2),
                 rank_from_top = rank(-x4$Chla)[x4$Unidad == "P6S1"], n_campaign = nrow(x4),
                 other_leaves_same_plot_Chla = paste(round(x4$Chla[x4$plot == 6 & x4$Unidad != "P6S1"], 3), collapse = " / "))
write_table(a1, "audit_leaves_per_campaign.csv"); write_table(a2, "audit_duplicated_keys.csv")
write_table(a3, "audit_missing_leaves.csv");       write_table(a4, "audit_block_by_treatment.csv")
write_table(a5, "audit_dropped_observation_P6S1_150DDC.csv")
cat("\nBlock x treatment (a balanced RCBD has a 1 in every cell):\n"); print(a4, row.names = FALSE)
cat("\nDropped by the original script (row 468):\n"); print(a5, row.names = FALSE)
